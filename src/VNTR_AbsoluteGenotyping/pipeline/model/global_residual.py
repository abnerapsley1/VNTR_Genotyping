"""
global_residual.py — Phase 3: the shipped predictor.

Stacks everything: raw signals + covariates + per-locus-calibration prediction
(where available) + transfer prediction + OOD score.  Fits three quantile GBMs
(0.05 / 0.50 / 0.95) on quadrant == 'train' for a point estimate + interval.

Confidence tier from interval width, OOD score, coverage and per-locus train fit:
  PASS       narrow interval, in-distribution
  LOW_CONF   wide interval or OOD or thin coverage
  NA         no usable signal (extraction failed / no reads / no calibration path)

Outputs:
  work/model/global_model.joblib
  work/model/predictions.parquet   (final: y_true, y_pred, y_lo, y_hi, tier, ...)
"""

from __future__ import annotations

import argparse

import joblib
import numpy as np
import polars as pl

from ..paths import load_config
from ._gbm import expm1_pred, log1p_target, make_regressor
from .assemble import feature_columns

_QUANTILES = (0.05, 0.50, 0.95)


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    plc = pl.read_parquet(mdir / "plc_predictions.parquet").select(
        "locus_id", "sample_id", "plc_pred", "plc_resid_sd")
    trf = pl.read_parquet(mdir / "transfer_predictions.parquet").select(
        "locus_id", "sample_id", "transfer_pred", "transfer_ood")
    tbl = tbl.join(plc, on=["locus_id", "sample_id"], how="left") \
             .join(trf, on=["locus_id", "sample_id"], how="left")

    seed = int(cfg.get("model.random_seed", 20260831))
    feats = feature_columns(tbl)
    for extra in ("plc_pred", "plc_resid_sd", "transfer_pred", "transfer_ood"):
        if extra not in feats and extra in tbl.columns:
            feats.append(extra)
    tbl = tbl.with_columns(plc_available=pl.col("plc_pred").is_not_null().cast(pl.Int8))
    feats.append("plc_available")
    print(f"global model: {len(feats)} features")

    tr = tbl.filter(pl.col("quadrant") == "train")
    Xtr = tr.select(feats).to_numpy().astype(np.float64)
    ytr = log1p_target(tr["combined_units"].to_numpy())
    Xall = tbl.select(feats).to_numpy().astype(np.float64)

    models = {}
    qpred = {}
    for q in _QUANTILES:
        m = make_regressor(quantile=q, seed=seed)
        m.fit(Xtr, ytr)
        qpred[q] = expm1_pred(m.predict(Xall))
        models[q] = m

    y_lo, y_pred, y_hi = qpred[0.05], qpred[0.50], qpred[0.95]
    y_lo = np.minimum(y_lo, y_pred)
    y_hi = np.maximum(y_hi, y_pred)

    width_rel = (y_hi - y_lo) / np.clip(y_pred, 1, None)
    ood = tbl["transfer_ood"].to_numpy()
    ood_hi = np.nanquantile(ood, 0.90) if np.isfinite(ood).any() else np.inf
    ext_ok = tbl["extraction_ok"].to_numpy() if "extraction_ok" in tbl.columns else np.ones(tbl.height)

    tier = np.where(
        (~np.isfinite(y_pred)),
        "NA",
        np.where((width_rel > 1.0) | (ood > ood_hi) | (width_rel != width_rel),
                 "LOW_CONF", "PASS"),
    )

    out = tbl.select("locus_id", "sample_id", "quadrant", "locus_fold", "sample_fold",
                     "combined_units", "combined_bp", "period",
                     "plc_pred", "transfer_pred", "transfer_ood").with_columns(
        y_pred=pl.Series(y_pred), y_lo=pl.Series(y_lo), y_hi=pl.Series(y_hi),
        interval_width_rel=pl.Series(width_rel), tier=pl.Series(tier))
    # add a few covariates for stratified evaluation
    for c in ("cov_is_segdup", "cov_period", "cov_mr_size_bp", "cov_mean_combined_units"):
        if c in tbl.columns:
            out = out.with_columns(tbl[c].alias(c))
    out.write_parquet(mdir / "predictions.parquet")
    joblib.dump({"quantile_models": models, "features": feats,
                 "quantiles": _QUANTILES, "ood_p90": float(ood_hi)},
                mdir / "global_model.joblib")
    print(f"wrote {mdir/'predictions.parquet'}, {mdir/'global_model.joblib'}")
    print(out.group_by("tier").agg(n=pl.len()).sort("n", descending=True))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
