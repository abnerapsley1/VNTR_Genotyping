"""
per_locus_calibration.py — Phase 1: one empirical calibration curve per locus.

For every locus with >= model.min_panel_samples_per_locus training-sample rows,
fit  combined_units = g(signal)  where g is isotonic (monotone) within the
observed signal range and robust-linear (Theil-Sen) for extrapolation.  The
1-D signal is chosen per locus as the highest train-correlation column from a
shortlist of normalised depth / IRR signals.

Generalisation axis here is SAMPLES (sample_fold): the locus is known, the
individual is new.  Novel-locus generalisation is the transfer model's job, so
only loci in the train locus-fold are calibrated (held-out loci get null).

Train-sample rows get OUT-OF-FOLD predictions (model.stack_cv_folds folds over
train samples): the global model stacks on plc_pred, and an in-sample curve
fits the panel almost exactly, which teaches the stacker to copy it.

Outputs:
  work/model/per_locus_calibration.parquet  (per-locus fit params + train stats)
  work/model/plc_predictions.parquet        (locus, sample, y_true, plc_pred, ...)
"""

from __future__ import annotations

import argparse
import hashlib
import warnings

import numpy as np
import polars as pl
from joblib import Parallel, delayed
from scipy.stats import spearmanr, theilslopes
from sklearn.isotonic import IsotonicRegression

from ..paths import load_config

try:                                   # the combined (ridge) curve lives in the vntrlen package
    from vntrlen import curves as _vc
except ImportError:                    # pragma: no cover
    _vc = None

_SIGNAL_SHORTLIST = [
    "depth_ratio_d20_w2000", "depth_ratio_d20_w1000", "depth_ratio_d20_w5000",
    "depth_ratio_d1_w2000", "array_depth_d20_norm",
    "irr_count_norm", "partial_repeat_frac", "n_spanning_norm",
    # added 2026-10-07: out-of-fold within-locus r 0.685 -> 0.714 on 30 samples
    "spanning_tlen_shift", "depth_ratio_d0_w2000", "anchored_irr_count_norm",
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


def _sample_cv_fold(sample_id: str, seed: int, k: int) -> int:
    h = hashlib.sha1(f"{seed}:cv:{sample_id}".encode()).digest()
    return int.from_bytes(h[:8], "big") % k


def _best_fit(x_cols: dict[str, np.ndarray], y: np.ndarray, min_n: int):
    """Pick the shortlist signal with the highest |train r|; return (col, params)."""
    if y.size < min_n:
        return None
    best = None
    for c, x in x_cols.items():
        params = _fit_one(x, y)
        if params is None:
            continue
        r = params["train_pearson"]
        score = abs(r) if np.isfinite(r) else -1
        if best is None or score > best[0]:
            best = (score, c, params)
    return None if best is None else best[1:]


def _fit(sig, y, min_n, method):
    if method == "ridge":
        return _vc.fit_ridge(sig, y, [c for c in _SIGNAL_SHORTLIST if c in sig], min_n)
    return _best_fit(sig, y, min_n)


def _pred(entry, sig, rows):
    if entry[0] == "ridge":
        return _vc.predict_entry(entry, {c: x[rows] for c, x in sig.items()})[0]
    return _predict_one(entry[1], sig[entry[0]][rows])


def _calibrate_locus(lid, sig: dict[str, np.ndarray], y: np.ndarray,
                     is_tr: np.ndarray, cvf: np.ndarray, k: int, min_n: int, method: str = "ridge"):
    """Fit one locus.  Returns (fit summary, preds, signal per row).

    Non-train-sample rows are predicted by the fit on all train samples; each
    train-sample row is predicted by a fit that excluded its own CV fold, so
    the global model never stacks on an in-sample calibration.
    """
    pred = np.full(y.size, np.nan)
    used = np.full(y.size, None, dtype=object)
    full = _fit({c: x[is_tr] for c, x in sig.items()}, y[is_tr], min_n, method)
    if full is None:
        return dict(locus_id=lid, signal=None, n=int(is_tr.sum()), train_pearson=None,
                    resid_sd=None, lin_slope=None, lin_intercept=None,
                    x_min=None, x_max=None), pred, used, np.nan
    col, params = full
    out = ~is_tr
    pred[out] = _pred(full, sig, out)
    used[out] = col
    for f in range(k):
        held = is_tr & (cvf == f)
        if not held.any():
            continue
        fit = is_tr & (cvf != f)
        sub = _fit({c: x[fit] for c, x in sig.items()}, y[fit], min_n, method)
        if sub is not None:
            pred[held] = _pred(sub, sig, held)
            used[held] = sub[0]
    summary = dict(locus_id=lid, signal=col, **{
        k_: params[k_] for k_ in ("n", "train_pearson", "resid_sd",
                                  "lin_slope", "lin_intercept", "x_min", "x_max")})
    return summary, pred, used, params["resid_sd"]


def _calibrate_chunk(groups, k, min_n, method="ridge"):
    warnings.filterwarnings("ignore")
    np.seterr(all="ignore")
    return [_calibrate_locus(*g, k=k, min_n=min_n, method=method) for g in groups]


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    min_n = int(cfg.get("model.min_panel_samples_per_locus", 8))
    k = int(cfg.get("model.stack_cv_folds", 5))
    seed = int(cfg.get("model.random_seed", 20260831))
    n_jobs = int(cfg.get("model.n_jobs", 16))
    method = str(cfg.get("model.curve_method", "ridge"))
    if method == "ridge" and _vc is None:
        raise SystemExit("model.curve_method 'ridge' needs the vntrlen package (pip install -e src/vntrlen)")

    cols_present = [c for c in _SIGNAL_SHORTLIST if c in tbl.columns]
    if not cols_present:
        raise SystemExit("none of the shortlist signal columns are in model_table")

    # Only loci in the train locus-fold get a panel calibration: a novel locus
    # (user BED, or held-out val/test loci) has no panel to calibrate against.
    keep = ["locus_id", "sample_id", "sample_fold", "locus_fold", "quadrant", "combined_units"]
    cal = tbl.filter(pl.col("locus_fold") == "train").select(keep + cols_present)
    cal = cal.with_columns(cvf=pl.col("sample_id").map_elements(
        lambda s: _sample_cv_fold(s, seed, k), return_dtype=pl.Int64))

    groups = []
    for (lid,), g in cal.group_by(["locus_id"], maintain_order=True):
        groups.append((lid, {c: g[c].to_numpy().astype(float) for c in cols_present},
                       g["combined_units"].to_numpy().astype(float),
                       (g["sample_fold"] == "train").to_numpy(),
                       g["cvf"].to_numpy()))
    n_chunks = max(1, n_jobs * 8)
    chunks = [groups[i::n_chunks] for i in range(n_chunks)]
    print(f"per-locus calibration ({method} curves): {len(groups):,} train-fold loci, {k}-fold cross-fit, "
          f"{n_jobs} jobs")
    results = [r for part in Parallel(n_jobs=n_jobs)(
        delayed(_calibrate_chunk)(c, k, min_n, method) for c in chunks) for r in part]

    by_lid = {g[0]: r for g, r in zip((g for c in chunks for g in c), results)}
    fits, preds = [], []
    for (lid,), g in cal.group_by(["locus_id"], maintain_order=True):
        summary, pred, used, rsd = by_lid[lid]
        fits.append(summary)
        preds.append(g.select(keep).with_columns(
            plc_pred=pl.Series(pred, dtype=pl.Float64),
            plc_signal=pl.Series([u for u in used], dtype=pl.Utf8),
            plc_resid_sd=pl.lit(rsd, dtype=pl.Float64)))

    fit_df = pl.DataFrame(fits, infer_schema_length=None)
    fit_df.write_parquet(mdir / "per_locus_calibration.parquet")
    # non-train-fold loci: present with null calibration (novel-locus path)
    rest = tbl.filter(pl.col("locus_fold") != "train").select(keep).with_columns(
        plc_pred=pl.lit(None, dtype=pl.Float64), plc_signal=pl.lit(None, dtype=pl.Utf8),
        plc_resid_sd=pl.lit(None, dtype=pl.Float64))
    pred_df = pl.concat(preds + [rest], how="diagonal_relaxed")
    pred_df.write_parquet(mdir / "plc_predictions.parquet")

    n_fit = fit_df.filter(pl.col("signal").is_not_null()).height
    print(f"per-locus calibration: fit {n_fit:,}/{fit_df.height:,} train-fold loci "
          f"(>= {min_n} train-sample rows)")
    _quick_metrics(pred_df)
    print(f"wrote {mdir/'per_locus_calibration.parquet'}, {mdir/'plc_predictions.parquet'}")


def _quick_metrics(pred: pl.DataFrame) -> None:
    d = pred.filter(pl.col("plc_pred").is_finite()
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
