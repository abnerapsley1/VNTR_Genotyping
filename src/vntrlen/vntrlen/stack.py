"""The cross-repeat model stack: transfer model + correction (point and 5%/95%) models
+ conformal interval shifts.  Applied after each repeat's calibration curve.

stack/ directory:
  stack.json        feature lists, conformal shifts per core-spread group, novelty-score data
  transfer.json     XGBoost booster: log1p(length) from reference-derivable features
  point.json        correction of log1p(anchor): squared error
  q_lo.json/q_hi.json   5th / 95th percentile of the correction
  ood.npz           standardised training-repeat covariates for the novelty score
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
from sklearn.neighbors import NearestNeighbors

from . import gbm

ALPHA = 0.10
CORE_SPREAD_EDGES = [1, 2, 5, 20, 100]
CORE_SPREAD_LABELS = ["0", "1", "2-4", "5-19", "20-99", ">=100"]


class Stack:
    def __init__(self, d: Path):
        self.dir = Path(d)
        self.meta = json.loads((self.dir / "stack.json").read_text())
        self.transfer = gbm.load(self.dir / "transfer.json")
        self.point = gbm.load(self.dir / "point.json")
        self.q_lo = gbm.load(self.dir / "q_lo.json")
        self.q_hi = gbm.load(self.dir / "q_hi.json")
        o = np.load(self.dir / "ood.npz")
        self.ood_mu, self.ood_sd, self.ood_T = o["mu"], o["sd"], o["T"]
        self._nn = None

    @property
    def features(self) -> list[str]:
        return self.meta["features"]

    @property
    def transfer_features(self) -> list[str]:
        return self.meta["transfer_features"]

    def ood_scores(self, cov: pl.DataFrame) -> np.ndarray:
        """Mean distance to the 5 nearest training repeats in standardised covariate space."""
        cols = self.meta["ood_columns"]
        if self._nn is None:
            self._nn = NearestNeighbors(n_neighbors=min(5, len(self.ood_T))).fit(self.ood_T)
        A = cov.select(cols).to_numpy().astype(float)
        A = np.nan_to_num((A - self.ood_mu) / self.ood_sd, nan=0.0)
        d, _ = self._nn.kneighbors(A)
        return d.mean(axis=1)


def _X(df: pl.DataFrame, cols: list[str]) -> np.ndarray:
    have = set(df.columns)
    return np.column_stack([df[c].cast(pl.Float64).to_numpy() if c in have else np.full(df.height, np.nan)
                            for c in cols])


def interval(lo, hi, t, lp):
    """Whole-unit interval: bounds moved out by t on the log1p scale (t < 0 narrows),
    floored/ceiled, and always containing the rounded point estimate (so a narrowing
    shift can never invert a narrow interval)."""
    pt = np.round(np.clip(np.expm1(lp), 0, None))
    a = np.floor(np.clip(np.expm1(lo - t), 0, None) + 1e-6)
    b = np.ceil(np.clip(np.expm1(hi + t), 0, None) - 1e-6)
    return np.minimum(a, pt), np.maximum(b, pt)


def core_groups(df: pl.DataFrame) -> np.ndarray:
    core = df["cov_panel_core90_units"].cast(pl.Float64).to_numpy() if "cov_panel_core90_units" in df.columns \
        else np.zeros(df.height)
    return np.array(CORE_SPREAD_LABELS)[np.digitize(np.nan_to_num(core, nan=0.0), CORE_SPREAD_EDGES)]


def add_transfer(df: pl.DataFrame, st: Stack) -> pl.DataFrame:
    tp = np.clip(np.expm1(gbm.predict(st.transfer, _X(df, st.transfer_features))), 0, None)
    return df.with_columns(transfer_pred=pl.Series(tp))


def anchor(df: pl.DataFrame) -> np.ndarray:
    a = df.select(pl.coalesce(pl.when(pl.col("plc_pred").is_finite()).then(pl.col("plc_pred")),
                              pl.col("cov_mean_combined_units"), 2 * pl.col("cov_ref_units_est"))
                  .clip(1, None))
    return a.to_series().to_numpy().astype(float)


def raw(df: pl.DataFrame, st: Stack):
    """(log1p point, log1p lower, log1p upper) before the conformal shift."""
    a = np.log1p(anchor(df))
    X = _X(df, st.features)
    p = gbm.predict(st.point, X)
    return a + p, a + np.minimum(gbm.predict(st.q_lo, X), p), a + np.maximum(gbm.predict(st.q_hi, X), p)


def shift_for(lo, hi, y, lp) -> float:
    """Smallest log1p shift whose whole-unit interval covers the conformal target."""
    n = y.size
    target = min(1.0, (1 - ALPHA) * (n + 1) / n)
    for t in np.round(np.arange(-0.5, 1.0001, 0.005), 3):
        a, b = interval(lo, hi, t, lp)
        if np.mean((y >= a) & (y <= b)) >= target:
            return float(t)
    return 1.0


def conformal_shifts(lo, hi, y, lp, groups, min_group: int = 200) -> tuple[float, dict]:
    t_all = shift_for(lo, hi, y, lp)
    return t_all, {g: (shift_for(lo[groups == g], hi[groups == g], y[groups == g], lp[groups == g])
                       if (groups == g).sum() >= min_group else t_all) for g in np.unique(groups)}


def predict(df: pl.DataFrame, st: Stack, t_all: float, shifts: dict) -> pl.DataFrame:
    lp, lo, hi = raw(df, st)
    g = core_groups(df)
    t = np.array([shifts.get(x, t_all) for x in g])
    y_lo, y_hi = interval(lo, hi, t, lp)
    return df.with_columns(est_units=pl.Series(np.clip(np.expm1(lp), 0, None)),
                           ci90_lo_units=pl.Series(y_lo), ci90_hi_units=pl.Series(y_hi),
                           interval_group=pl.Series(g))
