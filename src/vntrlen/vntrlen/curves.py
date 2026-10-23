"""Per-repeat calibration curves (short-read signals -> combined length).

Default ("ridge"): the shortlist signals are combined per repeat by ridge regression
(missing values filled with the training median, standardised, penalty chosen by
leave-one-out CV) into one score; the curve is then fitted on that score.
Alternative ("single"): the curve uses the one shortlist signal that best tracks the
truth (|Pearson r|).  Either way the curve is a monotone (isotonic) fit within the
training range and a robust Theil-Sen line beyond it.  Each repeat also stores a curve
built only from signals that exist for single-end reads.

A stored curve is (name, params): name is a signal (single) or "ridge".
"""
from __future__ import annotations

import json
import warnings

import numpy as np
from scipy.stats import theilslopes
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import RidgeCV

SIGNAL_SHORTLIST = [
    "depth_ratio_d20_w2000", "depth_ratio_d20_w1000", "depth_ratio_d20_w5000",
    "depth_ratio_d1_w2000", "array_depth_d20_norm",
    "irr_count_norm", "partial_repeat_frac", "n_spanning_norm",
    "spanning_tlen_shift", "depth_ratio_d0_w2000", "anchored_irr_count_norm",
]
PAIR_SIGNALS = {"n_spanning_norm", "spanning_tlen_shift", "anchored_irr_count_norm"}
SE_SHORTLIST = [s for s in SIGNAL_SHORTLIST if s not in PAIR_SIGNALS]
CURVE_METHOD = "ridge"
RIDGE_ALPHAS = np.logspace(-2, 3, 12)


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


def fit_ridge(sig: dict[str, np.ndarray], y: np.ndarray, shortlist: list[str], min_n: int):
    """("ridge", params): ridge-combined score of the shortlist signals + curve on that score."""
    ok = np.isfinite(y)
    if ok.sum() < min_n:
        return None
    names = [c for c in shortlist if c in sig]
    X = np.column_stack([sig[c] for c in names]).astype(float)[ok]
    yy = y[ok]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        med = np.nanmedian(X, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    Xf = np.where(np.isfinite(X), X, med)
    mu, sd = Xf.mean(0), Xf.std(0)
    sd = np.where(sd > 1e-9, sd, 1.0)
    base = dict(signals=names, med=med.tolist(), mu=mu.tolist(), sd=sd.tolist())
    if np.std(yy) == 0:                       # everyone has the same length at this repeat
        v = float(yy[0])
        return "ridge", dict(base, coef=[0.0] * len(names), intercept=v, x_min=v, x_max=v, iso_x=[v], iso_y=[v],
                             lin_slope=0.0, lin_intercept=v, resid_sd=0.0, n=int(yy.size), train_pearson=float("nan"))
    m = RidgeCV(alphas=RIDGE_ALPHAS).fit((Xf - mu) / sd, yy)
    p = dict(base, coef=m.coef_.tolist(), intercept=float(m.intercept_))
    score = _score(p, X)
    c = fit_one(score, yy)
    if c is None:
        return None
    return "ridge", dict(p, **c)


def _score(p: dict, X: np.ndarray) -> np.ndarray:
    X = np.where(np.isfinite(X), X, np.asarray(p["med"]))
    return ((X - np.asarray(p["mu"])) / np.asarray(p["sd"])) @ np.asarray(p["coef"]) + p["intercept"]


def predict_entry(entry, sig: dict[str, np.ndarray]):
    """(estimate, outside-training-range flag) for one stored curve on rows of signals."""
    name, p = entry
    if name == "ridge":
        n = len(next(iter(sig.values())))
        X = np.column_stack([sig.get(c, np.full(n, np.nan)) for c in p["signals"]]).astype(float)
        x = _score(p, X)
        if p["x_min"] == p["x_max"]:          # constant-length repeat
            return np.full(n, p["x_min"]), np.zeros(n, dtype=bool)
    else:
        x = np.asarray(sig[name], dtype=float)
    return predict_one(p, x), np.isfinite(x) & ((x < p["x_min"]) | (x > p["x_max"]))


def fit_entry(sig, y, shortlist, min_n, method: str | None = None):
    return (fit_ridge if (method or CURVE_METHOD) == "ridge" else best_fit)(sig, y, shortlist, min_n)


def fit_locus(sig: dict[str, np.ndarray], y: np.ndarray, min_n: int, method: str | None = None) -> dict:
    """Full-data curves for one repeat: {'pe': (name, params) | None, 'se': ...}."""
    return {"pe": fit_entry(sig, y, SIGNAL_SHORTLIST, min_n, method),
            "se": fit_entry(sig, y, SE_SHORTLIST, min_n, method)}


def oof_predict(sig: dict[str, np.ndarray], y: np.ndarray, folds: np.ndarray, min_n: int,
                single_end: np.ndarray | None = None, method: str | None = None):
    """Out-of-fold curve estimates: each row predicted by a curve fitted without its fold.
    Returns (pred, outside-range flag)."""
    pred = np.full(y.size, np.nan)
    outside = np.zeros(y.size, dtype=bool)
    se = np.zeros(y.size, dtype=bool) if single_end is None else single_end
    for f in np.unique(folds):
        held, fit = folds == f, folds != f
        sub = {c: v[fit] for c, v in sig.items()}
        for mode, rows in (("pe", held & ~se), ("se", held & se)):
            if not rows.any():
                continue
            e = fit_entry(sub, y[fit], SIGNAL_SHORTLIST if mode == "pe" else SE_SHORTLIST, min_n, method)
            if e is None:
                continue
            pred[rows], outside[rows] = predict_entry(e, {c: v[rows] for c, v in sig.items()})
    return pred, outside


def apply_curves(curves: dict, sig: dict[str, np.ndarray], single_end: bool):
    """Estimate for new people from a repeat's stored curves."""
    c = curves.get("se" if single_end else "pe") or curves.get("pe") or curves.get("se")
    n = len(next(iter(sig.values()))) if sig else 0
    if not c:
        return np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, dtype=bool), None
    pred, out = predict_entry(c, sig)
    return pred, np.full(n, c[1]["resid_sd"]), out, c[0]


def dumps(curves: dict) -> str:
    return json.dumps({k: (list(v) if v else None) for k, v in curves.items()})


def loads(s: str | None) -> dict:
    if not s:
        return {"pe": None, "se": None}
    d = json.loads(s)
    return {k: (tuple(v) if v else None) for k, v in d.items()}
