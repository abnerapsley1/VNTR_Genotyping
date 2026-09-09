"""
per_locus_calibration.py — Phase 1: one empirical calibration curve per locus.

For every locus with >= model.min_panel_samples_per_locus training-sample rows,
fit  combined_units = g(signal)  where g is isotonic (monotone) within the
observed signal range and robust-linear (Theil-Sen) for extrapolation.  The
1-D signal is chosen per locus as the highest train-correlation column from a
shortlist of normalised depth / IRR signals.

Generalisation axis here is SAMPLES (sample_fold): the locus is known, the
individual is new.  Novel-locus generalisation is the transfer model's job.

Outputs:
  work/model/per_locus_calibration.parquet  (per-locus fit params + train stats)
  work/model/plc_predictions.parquet        (locus, sample, y_true, plc_pred, ...)
"""

from __future__ import annotations

import argparse

import numpy as np
import polars as pl
from scipy.stats import spearmanr, theilslopes
from sklearn.isotonic import IsotonicRegression

from ..paths import load_config

_SIGNAL_SHORTLIST = [
    "depth_ratio_d20_w2000", "depth_ratio_d20_w1000", "depth_ratio_d20_w5000",
    "depth_ratio_d1_w2000", "array_depth_d20_norm",
    "irr_count_norm", "partial_repeat_frac", "n_spanning_norm",
]


def _fit_one(x: np.ndarray, y: np.ndarray):
    """Return dict of fit params or None if degenerate."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 4 or np.ptp(x) < 1e-9:
        return None
    order = np.argsort(x)
    x, y = x[order], y[order]
    iso = IsotonicRegression(increasing="auto", out_of_bounds="clip")
    iso.fit(x, y)
    yhat = iso.predict(x)
    resid = y - yhat
    try:
        slope, intercept, _, _ = theilslopes(y, x)
    except Exception:  # noqa: BLE001
        slope, intercept = np.polyfit(x, y, 1)
    r = float(np.corrcoef(x, y)[0, 1]) if x.size > 2 else np.nan
    return dict(
        x_min=float(x[0]), x_max=float(x[-1]),
        iso_x=iso.f_.x.tolist(), iso_y=iso.f_.y.tolist(),
        lin_slope=float(slope), lin_intercept=float(intercept),
        resid_sd=float(np.std(resid)), n=int(x.size),
        train_pearson=r,
    )


def _predict_one(params: dict, x: np.ndarray) -> np.ndarray:
    xi = np.asarray(params["iso_x"]); yi = np.asarray(params["iso_y"])
    out = np.interp(x, xi, yi, left=yi[0], right=yi[-1])
    lo = x < params["x_min"]
    hi = x > params["x_max"]
    lin = params["lin_slope"] * x + params["lin_intercept"]
    out = np.where(lo | hi, lin, out)
    return np.clip(out, 0, None)


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    min_n = int(cfg.get("model.min_panel_samples_per_locus", 8))

    cols_present = [c for c in _SIGNAL_SHORTLIST if c in tbl.columns]
    if not cols_present:
        raise SystemExit("none of the shortlist signal columns are in model_table")

    fits: list[dict] = []
    preds: list[pl.DataFrame] = []
    tr_mask = pl.col("sample_fold") == "train"

    for (lid,), g in tbl.group_by(["locus_id"], maintain_order=True):
        gtr = g.filter(tr_mask)
        y_tr = gtr["combined_units"].to_numpy()
        best = None
        if gtr.height >= min_n:
            for c in cols_present:
                x_tr = gtr[c].to_numpy()
                params = _fit_one(x_tr, y_tr)
                if params is None:
                    continue
                score = abs(params["train_pearson"]) if np.isfinite(params["train_pearson"]) else -1
                if best is None or score > best[0]:
                    best = (score, c, params)
        if best is None:
            fits.append(dict(locus_id=lid, signal=None, n=gtr.height, train_pearson=None,
                             resid_sd=None, lin_slope=None, lin_intercept=None,
                             x_min=None, x_max=None))
            preds.append(g.select("locus_id", "sample_id", "sample_fold", "locus_fold",
                                  "quadrant", "combined_units").with_columns(
                plc_pred=pl.lit(None, dtype=pl.Float64),
                plc_signal=pl.lit(None, dtype=pl.Utf8),
                plc_resid_sd=pl.lit(None, dtype=pl.Float64)))
            continue

        _, sigcol, params = best
        x_all = g[sigcol].to_numpy()
        yhat = _predict_one(params, x_all)
        fits.append(dict(locus_id=lid, signal=sigcol, **{
            k: params[k] for k in ("n", "train_pearson", "resid_sd",
                                   "lin_slope", "lin_intercept", "x_min", "x_max")}))
        preds.append(g.select("locus_id", "sample_id", "sample_fold", "locus_fold",
                              "quadrant", "combined_units").with_columns(
            plc_pred=pl.Series(yhat),
            plc_signal=pl.lit(sigcol),
            plc_resid_sd=pl.lit(params["resid_sd"])))

    fit_df = pl.DataFrame(fits)
    fit_df.write_parquet(mdir / "per_locus_calibration.parquet")
    pred_df = pl.concat(preds, how="diagonal_relaxed")
    pred_df.write_parquet(mdir / "plc_predictions.parquet")

    n_fit = fit_df.filter(pl.col("signal").is_not_null()).height
    print(f"per-locus calibration: fit {n_fit:,}/{fit_df.height:,} loci "
          f"(>= {min_n} train-sample rows)")
    _quick_metrics(pred_df)
    print(f"wrote {mdir/'per_locus_calibration.parquet'}, {mdir/'plc_predictions.parquet'}")


def _quick_metrics(pred: pl.DataFrame) -> None:
    d = pred.filter(pl.col("plc_pred").is_not_null()
                    & pl.col("combined_units").is_not_null())
    for q in ["train", "test_novel_sample", "test_novel_locus", "test_novel_both", "val"]:
        s = d.filter(pl.col("quadrant") == q)
        if s.height < 5:
            continue
        yt = s["combined_units"].to_numpy()
        yp = s["plc_pred"].to_numpy()
        mae = float(np.median(np.abs(yp - yt)))
        rho = spearmanr(yp, yt).statistic
        rel = float(np.median(np.abs(yp - yt) / np.clip(yt, 1, None)))
        print(f"  {q:18s} n={s.height:>7,}  MAE={mae:6.2f}u  medRelErr={rel:5.2f}  "
              f"spearman={rho:+.3f}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
