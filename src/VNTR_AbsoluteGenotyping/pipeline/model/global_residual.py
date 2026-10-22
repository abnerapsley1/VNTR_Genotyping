"""
global_residual.py — Phase 3: the shipped predictor.

Two paths, chosen per locus:
  known  bundled loci with a panel (locus_fold == train): raw signals + all
         covariates + cross-fitted per-locus-calibration and transfer
         predictions + OOD score.  Anchor = plc_pred, else the panel mean.
  novel  any other locus (held-out loci here; a user-BED locus in deployment):
         novel-safe features only — no panel statistics, no per-locus fit.
         Anchor = 2 x reference copy number.

Each path models the RESIDUAL  log1p(y) - log1p(anchor)  so the GBM spends its
capacity on within-locus (between-individual) variation instead of the
between-locus length scale (2 .. >1000 units).  Per path:
  point     squared-error GBM            -> y_pred
  interval  quantile GBMs (0.05 / 0.95)  -> raw interval, reported in whole
            motif units (floor/ceil) and shifted by a split-conformal step on
            the validation rows of that path, separately per group (known:
            locus core-spread bin; novel: predicted-width quintile), so the
            rounded interval covers 90% within groups, not only on average.

Confidence tier:
  PASS       calibrated interval narrow, in-distribution
  LOW_CONF   wide interval or OOD
  NA         no usable prediction

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
from ._gbm import make_regressor
from .assemble import feature_columns
from .calibration_transfer import novel_safe_features

_ALPHA = 0.10
_ANCHOR_COLS = ("anchor", "anchor_known", "anchor_novel")


def _X(df: pl.DataFrame, cols: list[str]) -> np.ndarray:
    return df.select(cols).to_numpy().astype(np.float64)


CORE_SPREAD_EDGES = [1, 2, 5, 20, 100]          # bins: 0 | 1 | 2-4 | 5-19 | 20-99 | >=100
CORE_SPREAD_LABELS = ["0", "1", "2-4", "5-19", "20-99", ">=100"]
_MIN_GROUP_CAL = 200


def _shift_for(lo_raw, hi_raw, y, interval) -> float:
    """Smallest log1p shift t whose whole-unit interval covers the conformal
    target on these calibration rows (coverage is monotone in t)."""
    n = y.size
    target = min(1.0, (1 - _ALPHA) * (n + 1) / n)
    for t in np.round(np.arange(-0.5, 1.0001, 0.005), 3):
        lo_c, hi_c = interval(lo_raw, hi_raw, t)
        if np.mean((y >= lo_c) & (y <= hi_c)) >= target:
            return float(t)
    return 1.0


def _fit_path(name: str, tr: pl.DataFrame, cal: pl.DataFrame, pred: pl.DataFrame,
              feats: list[str], anchor: str, group_mode: str,
              seed: int) -> tuple[dict, dict[str, np.ndarray]]:
    """Fit point + interval models on residual vs this path's `anchor` column
    (train, calibration and prediction rows all use the same anchor); then
    conformalize per group on `cal` (Mondrian split-conformal).

    group_mode "core":  locus core spread bins (needs cov_panel_core90_units;
                        known path only — a novel locus has no panel)
    group_mode "width": quintiles of the raw predicted interval width
    """
    def resid(df):
        return np.log1p(df["combined_units"].to_numpy()) - np.log1p(df[anchor].to_numpy())

    r_tr = resid(tr)
    point = make_regressor(seed=seed).fit(_X(tr, feats), r_tr)
    lo_m = make_regressor(quantile=_ALPHA / 2, seed=seed).fit(_X(tr, feats), r_tr)
    hi_m = make_regressor(quantile=1 - _ALPHA / 2, seed=seed).fit(_X(tr, feats), r_tr)

    def raw(df):
        """(log1p point, log1p lower, log1p upper) before any conformal shift."""
        a = np.log1p(df[anchor].to_numpy())
        X_ = _X(df, feats)
        p_ = point.predict(X_)
        return a + p_, a + np.minimum(lo_m.predict(X_), p_), a + np.maximum(hi_m.predict(X_), p_)

    def interval(lo_, hi_, t):
        """Whole-unit interval: bounds moved out by t (log1p), floored/ceiled."""
        return (np.floor(np.clip(np.expm1(lo_ - t), 0, None) + 1e-6),
                np.ceil(np.clip(np.expm1(hi_ + t), 0, None) - 1e-6))

    width_edges: list[float] = []

    def groups(df, lo_, hi_) -> np.ndarray:
        if group_mode == "core":
            core = df["cov_panel_core90_units"].to_numpy().astype(float)
            g = np.digitize(np.nan_to_num(core, nan=0.0), CORE_SPREAD_EDGES)
            return np.array(CORE_SPREAD_LABELS)[g]
        g = np.digitize(hi_ - lo_, width_edges)
        return np.array([f"width_q{k + 1}" for k in range(len(width_edges) + 1)])[g]

    # Conformal step in whole units, per group: lengths are integer motif
    # counts, so the interval is reported as integers; each group gets the
    # smallest shift whose rounded interval covers the target on its held-out
    # calibration rows, so coverage holds within groups, not just on average.
    shifts: dict[str, float] = {}
    t_all = 0.0
    if cal.height >= 50:
        y_cal = cal["combined_units"].to_numpy()
        _, lo_c, hi_c = raw(cal)
        if group_mode == "width":
            width_edges = list(np.quantile(hi_c - lo_c, [0.2, 0.4, 0.6, 0.8]))
        t_all = _shift_for(lo_c, hi_c, y_cal, interval)
        g_cal = groups(cal, lo_c, hi_c)
        for g in np.unique(g_cal):
            m = g_cal == g
            shifts[g] = (_shift_for(lo_c[m], hi_c[m], y_cal[m], interval)
                         if m.sum() >= _MIN_GROUP_CAL else t_all)
    print(f"  {name}: train {tr.height:,} rows, conformal calibration on {cal.height:,} "
          f"rows; shift overall {t_all:+.3f}, per group "
          + ", ".join(f"{k}:{v:+.3f}" for k, v in shifts.items()))

    lp, lo_raw, hi_raw = raw(pred)
    g_pred = groups(pred, lo_raw, hi_raw)
    t_pred = np.array([shifts.get(g, t_all) for g in g_pred])
    y_lo, y_hi = interval(lo_raw, hi_raw, t_pred)
    out = {"y_pred": np.clip(np.expm1(lp), 0, None), "y_lo": y_lo, "y_hi": y_hi,
           "interval_group": g_pred}
    models = {"point": point, "q_lo": lo_m, "q_hi": hi_m, "conformal_shift": t_all,
              "conformal_shift_by_group": shifts, "group_mode": group_mode,
              "width_edges": width_edges, "features": feats}
    return models, out


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    plc = pl.read_parquet(mdir / "plc_predictions.parquet").select(
        "locus_id", "sample_id", "plc_pred", "plc_resid_sd")
    trf = pl.read_parquet(mdir / "transfer_predictions.parquet").select(
        "locus_id", "sample_id", "transfer_pred", "transfer_ood")
    tbl = tbl.join(plc, on=["locus_id", "sample_id"], how="left") \
             .join(trf, on=["locus_id", "sample_id"], how="left")
    tbl = tbl.with_columns(
        plc_pred=pl.when(pl.col("plc_pred").is_finite()).then(pl.col("plc_pred")))
    tbl = tbl.with_columns(plc_available=pl.col("plc_pred").is_not_null().cast(pl.Int8))
    seed = int(cfg.get("model.random_seed", 20260831))

    is_known = pl.col("locus_fold") == "train"
    ref2 = 2 * pl.col("cov_ref_units_est")
    tbl = tbl.with_columns(
        anchor_known=pl.coalesce("plc_pred", "cov_mean_combined_units", ref2).clip(1, None),
        anchor_novel=ref2.clip(1, None))

    known_feats = [c for c in dict.fromkeys(feature_columns(tbl)) if c not in _ANCHOR_COLS]
    novel_feats = [c for c in novel_safe_features(tbl)
                   if c not in _ANCHOR_COLS and c != "transfer_pred"
                   and not c.startswith("transfer_")]
    print(f"global model: known path {len(known_feats)} features, "
          f"novel path {len(novel_feats)} features")

    tr = tbl.filter(pl.col("quadrant") == "train")
    # conformal calibration rows: validation rows that look like each path's use
    cal_known = tbl.filter(is_known & (pl.col("sample_fold") == "val"))
    cal_novel = tbl.filter((pl.col("locus_fold") == "val") & (pl.col("sample_fold") == "train"))
    known_rows = tbl.filter(is_known)
    novel_rows = tbl.filter(~is_known)

    km, kp = _fit_path("known", tr, cal_known, known_rows, known_feats, "anchor_known",
                       "core", seed)
    nm, npred = _fit_path("novel", tr, cal_novel, novel_rows, novel_feats, "anchor_novel",
                          "width", seed)

    parts = []
    for rows, p, path in ((known_rows, kp, "known"), (novel_rows, npred, "novel")):
        parts.append(rows.select(
            "locus_id", "sample_id", "quadrant", "locus_fold", "sample_fold",
            "combined_units", "combined_bp", "period", pl.col(f"anchor_{path}").alias("anchor"),
            "plc_pred", "transfer_pred", "transfer_ood",
            *[c for c in ("cov_is_segdup", "cov_period", "cov_mr_size_bp",
                          "cov_mean_combined_units", "cov_panel_median_units",
                          "cov_panel_p5_units", "cov_panel_p95_units",
                          "cov_panel_core90_units", "cov_panel_n_outliers")
              if c in rows.columns],
        ).with_columns(y_pred=pl.Series(p["y_pred"]), y_lo=pl.Series(p["y_lo"]),
                       y_hi=pl.Series(p["y_hi"]), path=pl.lit(path),
                       interval_group=pl.Series(p["interval_group"])))
    out = pl.concat(parts)

    width_max = float(cfg.get("model.pass_max_interval_width_rel", 1.0))
    ood = out["transfer_ood"].to_numpy()
    ood_hi = np.nanquantile(ood, 0.90) if np.isfinite(ood).any() else np.inf
    y_pred, y_lo, y_hi = (out[c].to_numpy() for c in ("y_pred", "y_lo", "y_hi"))
    width_rel = (y_hi - y_lo) / np.clip(y_pred, 1, None)
    tier = np.where(
        ~np.isfinite(y_pred), "NA",
        np.where((width_rel > width_max) | (ood > ood_hi) | ~np.isfinite(width_rel),
                 "LOW_CONF", "PASS"))
    out = out.with_columns(interval_width_rel=pl.Series(width_rel), tier=pl.Series(tier))

    out.write_parquet(mdir / "predictions.parquet")
    joblib.dump({"known": km, "novel": nm, "alpha": _ALPHA, "ood_p90": float(ood_hi),
                 "pass_max_interval_width_rel": width_max},
                mdir / "global_model.joblib")
    print(f"wrote {mdir/'predictions.parquet'}, {mdir/'global_model.joblib'}")
    print(out.group_by("path", "tier").agg(n=pl.len()).sort("path", "tier"))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
