"""
calibration_transfer.py — Phase 2: predict combined length for ANY locus from
signals + locus covariates, WITHOUT using locus identity.

This is the model that carries a novel user-BED locus (no panel calibration).
Trained on quadrant == 'train'; the headline metric is quadrant
'test_novel_locus' (held-out loci, known samples).

Also fits a simple out-of-distribution score: a locus far (in covariate space)
from the training loci gets its confidence downgraded downstream.

Outputs:
  work/model/transfer_model.joblib
  work/model/transfer_predictions.parquet   (locus, sample, transfer_pred, ood_score)
"""

from __future__ import annotations

import argparse

import joblib
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

from ..paths import load_config
from ._gbm import expm1_pred, log1p_target, make_regressor
from .assemble import feature_columns

# columns that would leak locus identity / panel truth into the transfer model
_LEAK = {"cov_mean_combined_units", "cov_var_combined_units", "cov_min_combined_units",
         "cov_max_combined_units", "cov_mean_combined_bp", "cov_var_combined_bp",
         "cov_panel_cn_range_units", "cov_n_called", "cov_h_exp", "cov_h_obs",
         "cov_n_distinct_hap_units", "cov_reciprocal_overlap", "cov_matched",
         "cov_ref_units_est"}


def _matrix(tbl: pl.DataFrame, cols: list[str]) -> np.ndarray:
    return tbl.select(cols).to_numpy().astype(np.float64)


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    seed = int(cfg.get("model.random_seed", 20260831))

    feats = [c for c in feature_columns(tbl) if c not in _LEAK]
    feats = [c for c in feats if not c.startswith("plc_")]
    print(f"transfer model: {len(feats)} features, target=log1p(combined_units)")

    tr = tbl.filter(pl.col("quadrant") == "train")
    if tr.height < 200:
        raise SystemExit(f"too few training rows ({tr.height}) — need more samples/loci")

    Xtr = _matrix(tr, feats)
    ytr = log1p_target(tr["combined_units"].to_numpy())
    model = make_regressor(seed=seed)
    model.fit(Xtr, ytr)

    Xall = _matrix(tbl, feats)
    pred = expm1_pred(model.predict(Xall))

    # OOD score: distance to nearest training locus in standardised covariate space
    cov_cols = [c for c in feats if c.startswith("cov_")]
    ood = _ood_scores(tbl, tr, cov_cols)

    out = tbl.select("locus_id", "sample_id", "quadrant", "locus_fold",
                     "sample_fold", "combined_units").with_columns(
        transfer_pred=pl.Series(pred),
        transfer_ood=pl.Series(ood))
    out.write_parquet(mdir / "transfer_predictions.parquet")
    joblib.dump({"model": model, "features": feats, "leak_excluded": sorted(_LEAK)},
                mdir / "transfer_model.joblib")

    _metrics(out)
    print(f"wrote {mdir/'transfer_model.joblib'}, {mdir/'transfer_predictions.parquet'}")


def _ood_scores(tbl, tr, cov_cols) -> np.ndarray:
    if not cov_cols:
        return np.zeros(tbl.height)
    loci_all = tbl.select(["locus_id"] + cov_cols).unique(subset="locus_id")
    loci_tr = tr.select(["locus_id"] + cov_cols).unique(subset="locus_id")
    A = loci_all.select(cov_cols).to_numpy().astype(float)
    T = loci_tr.select(cov_cols).to_numpy().astype(float)
    mu = np.nanmean(T, axis=0)
    sd = np.nanstd(T, axis=0) + 1e-9
    A = np.nan_to_num((A - mu) / sd, nan=0.0)
    T = np.nan_to_num((T - mu) / sd, nan=0.0)
    nn = NearestNeighbors(n_neighbors=min(5, len(T))).fit(T)
    dist, _ = nn.kneighbors(A)
    score = dist.mean(axis=1)
    per_locus = dict(zip(loci_all["locus_id"].to_list(), score))
    return np.array([per_locus.get(l, np.nan) for l in tbl["locus_id"].to_list()])


def _metrics(out: pl.DataFrame) -> None:
    for q in ["train", "test_novel_sample", "test_novel_locus", "test_novel_both", "val"]:
        s = out.filter((pl.col("quadrant") == q)
                       & pl.col("transfer_pred").is_not_null()
                       & pl.col("combined_units").is_not_null())
        if s.height < 5:
            continue
        yt = s["combined_units"].to_numpy()
        yp = s["transfer_pred"].to_numpy()
        mae = float(np.median(np.abs(yp - yt)))
        rel = float(np.median(np.abs(yp - yt) / np.clip(yt, 1, None)))
        rho = spearmanr(yp, yt).statistic
        # per-locus median spearman
        pl_rho = []
        for (_l,), gg in s.group_by(["locus_id"]):
            if gg.height < 6:
                continue
            a = gg["transfer_pred"].to_numpy()
            b = gg["combined_units"].to_numpy()
            if np.std(a) < 1e-9 or np.std(b) < 1e-9:
                continue
            pl_rho.append(spearmanr(a, b).statistic)
        mpl = float(np.nanmedian(pl_rho)) if pl_rho else float("nan")
        print(f"  {q:18s} n={s.height:>7,}  MAE={mae:6.2f}u  medRelErr={rel:5.2f}  "
              f"spearman(pool)={rho:+.3f}  median per-locus spearman={mpl:+.3f}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
