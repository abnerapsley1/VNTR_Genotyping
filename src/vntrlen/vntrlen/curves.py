"""Per-repeat calibration curves (one short-read signal -> combined length).

Within the training range: monotone (isotonic) fit; beyond it: robust Theil-Sen line.
The signal is the shortlist entry that best tracks the truth at that repeat.  Each
repeat also stores a curve restricted to signals that exist for single-end reads.
"""
from __future__ import annotations

import json
import warnings

import numpy as np
from scipy.stats import theilslopes
from sklearn.isotonic import IsotonicRegression

SIGNAL_SHORTLIST = [
    "depth_ratio_d20_w2000", "depth_ratio_d20_w1000", "depth_ratio_d20_w5000",
    "depth_ratio_d1_w2000", "array_depth_d20_norm",
    "irr_count_norm", "partial_repeat_frac", "n_spanning_norm",
    "spanning_tlen_shift", "depth_ratio_d0_w2000", "anchored_irr_count_norm",
]
PAIR_SIGNALS = {"n_spanning_norm", "spanning_tlen_shift", "anchored_irr_count_norm"}
SE_SHORTLIST = [s for s in SIGNAL_SHORTLIST if s not in PAIR_SIGNALS]


def fit_one(x: np.ndarray, y: np.ndarray) -> dict | None:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 4 or np.ptp(x) < 1e-9:
        return None
    o = np.argsort(x)
    x, y = x[o], y[o]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        iso = IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(x, y)
        try:
            slope, intercept, _, _ = theilslopes(y, x)
        except Exception:  # noqa: BLE001
            slope, intercept = np.polyfit(x, y, 1)
    resid = y - iso.predict(x)
    r = float(np.corrcoef(x, y)[0, 1]) if x.size > 2 and np.std(y) > 0 else np.nan
    return dict(x_min=float(x[0]), x_max=float(x[-1]), iso_x=iso.f_.x.tolist(), iso_y=iso.f_.y.tolist(),
                lin_slope=float(slope), lin_intercept=float(intercept), resid_sd=float(np.std(resid)),
                n=int(x.size), train_pearson=r)


def predict_one(p: dict, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    xi, yi = np.asarray(p["iso_x"]), np.asarray(p["iso_y"])
    out = np.interp(x, xi, yi, left=yi[0], right=yi[-1])
    out = np.where((x < p["x_min"]) | (x > p["x_max"]), p["lin_slope"] * x + p["lin_intercept"], out)
    return np.where(np.isfinite(x), np.clip(out, 0, None), np.nan)


def best_fit(sig: dict[str, np.ndarray], y: np.ndarray, shortlist: list[str], min_n: int):
    """(signal, params) for the shortlist signal with the highest |train r|, or None."""
    if np.isfinite(y).sum() < min_n:
        return None
    best = None
    for c in shortlist:
        if c not in sig:
            continue
        p = fit_one(sig[c], y)
        if p is None:
            continue
        score = abs(p["train_pearson"]) if np.isfinite(p["train_pearson"]) else -1
        if best is None or score > best[0]:
            best = (score, c, p)
    return None if best is None else (best[1], best[2])


def fit_locus(sig: dict[str, np.ndarray], y: np.ndarray, min_n: int) -> dict:
    """Full-data curves for one repeat: {'pe': (signal, params) | None, 'se': ...}."""
    return {"pe": best_fit(sig, y, SIGNAL_SHORTLIST, min_n), "se": best_fit(sig, y, SE_SHORTLIST, min_n)}


def oof_predict(sig: dict[str, np.ndarray], y: np.ndarray, folds: np.ndarray, min_n: int,
                single_end: np.ndarray | None = None):
    """Out-of-fold curve estimates: each row predicted by a curve fitted without its fold.
    Returns (pred, resid_sd of the full fit, signal used per row, outside-range flag)."""
    pred = np.full(y.size, np.nan)
    outside = np.zeros(y.size, dtype=bool)
    se = np.zeros(y.size, dtype=bool) if single_end is None else single_end
    for f in np.unique(folds):
        held, fit = folds == f, folds != f
        sub = {c: v[fit] for c, v in sig.items()}
        for mode, rows in (("pe", held & ~se), ("se", held & se)):
            if not rows.any():
                continue
            bf = best_fit(sub, y[fit], SIGNAL_SHORTLIST if mode == "pe" else SE_SHORTLIST, min_n)
            if bf is None:
                continue
            c, p = bf
            x = sig[c][rows]
            pred[rows] = predict_one(p, x)
            outside[rows] = (x < p["x_min"]) | (x > p["x_max"])
    return pred, outside


def apply_curves(curves: dict, sig: dict[str, np.ndarray], single_end: bool):
    """Estimate for new people from a repeat's stored curves."""
    c = curves.get("se" if single_end else "pe") or curves.get("pe") or curves.get("se")
    n = len(next(iter(sig.values()))) if sig else 0
    if not c:
        return np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, dtype=bool), None
    name, p = c
    x = sig[name]
    return predict_one(p, x), np.full(n, p["resid_sd"]), (x < p["x_min"]) | (x > p["x_max"]), name


def dumps(curves: dict) -> str:
    return json.dumps({k: (list(v) if v else None) for k, v in curves.items()})


def loads(s: str | None) -> dict:
    if not s:
        return {"pe": None, "se": None}
    d = json.loads(s)
    return {k: (tuple(v) if v else None) for k, v in d.items()}
