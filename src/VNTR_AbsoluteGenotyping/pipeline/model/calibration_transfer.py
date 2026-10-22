"""
calibration_transfer.py — Phase 2: predict combined length for ANY locus from
signals + locus covariates, WITHOUT using locus identity.

This is the model that carries a novel user-BED locus (no panel calibration).
It sees only NOVEL_SAFE_COVARIATES (BED/reference-derivable), never population
long-read statistics.  Trained on quadrant == 'train'; the headline metric is
quadrant 'test_novel_locus' (held-out loci, known samples).  Train-row
predictions are cross-fitted over samples for stacking.

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
from .per_locus_calibration import _sample_cv_fold

# Covariates a NOVEL locus (arbitrary user BED, not in the vamos catalogue) can
# actually have: derivable from the BED + reference sequence/annotation tracks.
# Everything else in locus_covariates is computed from population long-read
# allele data (panel length mean/var/range, heterozygosity, allele counts,
# anova population means, mutRate length/mutation stats) and would leak the
# answer into the novel-locus path.
NOVEL_SAFE_COVARIATES = {
    "cov_period", "cov_array_span_bp", "cov_log_array_span_bp", "cov_ref_units_est",
    "cov_consensus_gc", "cov_consensus_base_entropy", "cov_consensus_len",
    "cov_is_segdup", "cov_mr_segdup_identity", "cov_mr_segdup_cnv",
    "cov_mr_cons_gc", "cov_mr_all_gc", "cov_mr_te_frac_locus", "cov_mr_te_frac_element",
    "cov_mr_coding", "cov_mr_exon", "cov_mr_transcribed",
    "cov_mr_dist_to_recomb", "cov_mr_telomere", "cov_mr_size_bp", "cov_mr_period",
}


def novel_safe_features(tbl: pl.DataFrame) -> list[str]:
    """Read signals + nuis_* + the reference-derivable covariate allowlist."""
    return [c for c in feature_columns(tbl)
            if not c.startswith("plc_")
            and (not c.startswith("cov_") or c in NOVEL_SAFE_COVARIATES)]


def _matrix(tbl: pl.DataFrame, cols: list[str]) -> np.ndarray:
    return tbl.select(cols).to_numpy().astype(np.float64)


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    seed = int(cfg.get("model.random_seed", 20260831))

    feats = novel_safe_features(tbl)
    print(f"transfer model: {len(feats)} features (novel-safe), target=log1p(combined_units)")

    tr = tbl.filter(pl.col("quadrant") == "train")
    if tr.height < 200:
        raise SystemExit(f"too few training rows ({tr.height}) — need more samples/loci")

    model = make_regressor(seed=seed)
    model.fit(_matrix(tr, feats), log1p_target(tr["combined_units"].to_numpy()))
    pred = expm1_pred(model.predict(_matrix(tbl, feats)))

    # train rows: out-of-fold over train samples, so the stacked global model
    # sees transfer_pred as it behaves on an unseen individual
    k = int(cfg.get("model.stack_cv_folds", 5))
    cvf = np.array([_sample_cv_fold(s, seed, k) for s in tr["sample_id"].to_list()])
    tr_idx = np.flatnonzero((tbl["quadrant"] == "train").to_numpy())
    for f in range(k):
        held = cvf == f
        if not held.any() or held.all():
            continue
        m = make_regressor(seed=seed)
        m.fit(_matrix(tr.filter(pl.Series(~held)), feats),
              log1p_target(tr["combined_units"].to_numpy()[~held]))
        pred[tr_idx[held]] = expm1_pred(m.predict(_matrix(tr.filter(pl.Series(held)), feats)))
        print(f"  cross-fit fold {f + 1}/{k} done")

    # OOD score: distance to nearest training locus in standardised covariate space
    cov_cols = [c for c in feats if c.startswith("cov_")]
    ood = _ood_scores(tbl, tr, cov_cols)

    out = tbl.select("locus_id", "sample_id", "quadrant", "locus_fold",
                     "sample_fold", "combined_units").with_columns(
        transfer_pred=pl.Series(pred),
        transfer_ood=pl.Series(ood))
    out.write_parquet(mdir / "transfer_predictions.parquet")
    joblib.dump({"model": model, "features": feats,
                 "novel_safe_covariates": sorted(NOVEL_SAFE_COVARIATES)},
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
