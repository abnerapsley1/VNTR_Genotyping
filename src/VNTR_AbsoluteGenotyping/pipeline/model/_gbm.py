"""
_gbm.py — gradient-boosted regressor wrapper.

Uses sklearn HistGradientBoostingRegressor: native NaN handling, fast on ~1-10M
rows, supports quantile loss for prediction intervals.  If xgboost or lightgbm
is installed (GPU-capable), set VNTR_GBM=xgboost|lightgbm to use it instead.
"""

from __future__ import annotations

import os

import numpy as np


def _backend() -> str:
    return os.environ.get("VNTR_GBM", "sklearn").lower()


def make_regressor(*, quantile: float | None = None, seed: int = 0,
                   n_estimators: int = 600, learning_rate: float = 0.05,
                   max_depth: int | None = None, max_leaves: int = 63):
    b = _backend()
    if b == "xgboost":
        # GPU stand-in for the sklearn model below: same leaf-wise growth, leaf
        # cap, L2 penalty, bins and early-stopping rule, so only the backend differs
        import xgboost as xgb
        obj = "reg:quantileerror" if quantile is not None else "reg:squarederror"
        kw = dict(objective=obj, n_estimators=n_estimators, learning_rate=learning_rate,
                  grow_policy="lossguide", max_leaves=max_leaves, max_depth=max_depth or 0,
                  reg_lambda=1.0, max_bin=255, tree_method="hist",
                  device=os.environ.get("VNTR_GBM_DEVICE", "cuda"),
                  early_stopping_rounds=30, random_state=seed)
        if quantile is not None:
            kw["quantile_alpha"] = quantile
        return _XGBEarlyStop(xgb.XGBRegressor(**kw), seed)
    if b == "lightgbm":
        import lightgbm as lgb
        return lgb.LGBMRegressor(
            objective="quantile" if quantile is not None else "regression",
            alpha=quantile if quantile is not None else 0.9,
            n_estimators=n_estimators, learning_rate=learning_rate,
            num_leaves=max_leaves, subsample=0.8, subsample_freq=1,
            colsample_bytree=0.8, random_state=seed, n_jobs=-1)
    from sklearn.ensemble import HistGradientBoostingRegressor
    return HistGradientBoostingRegressor(
        loss="quantile" if quantile is not None else "squared_error",
        quantile=quantile,
        max_iter=n_estimators, learning_rate=learning_rate,
        max_depth=max_depth, max_leaf_nodes=max_leaves,
        l2_regularization=1.0, early_stopping=True, validation_fraction=0.1,
        n_iter_no_change=30, random_state=seed)


class _XGBEarlyStop:
    """XGBRegressor fit with sklearn's early-stopping rule: a random 10% of the
    training rows is held out and boosting stops after 30 rounds without gain."""

    def __init__(self, model, seed: int):
        self.model, self.seed = model, seed

    def fit(self, X, y):
        idx = np.random.default_rng(self.seed).permutation(len(y))
        n_val = max(1, int(0.1 * len(y)))
        va, tr = idx[:n_val], idx[n_val:]
        self.model.fit(X[tr], y[tr], eval_set=[(X[va], y[va])], verbose=False)
        return self

    def predict(self, X):
        return self.model.predict(X)


def fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xall: np.ndarray, **kw) -> np.ndarray:
    m = make_regressor(**kw)
    m.fit(Xtr, ytr)
    return m.predict(Xall)


def log1p_target(y: np.ndarray) -> np.ndarray:
    return np.log1p(np.clip(y, 0, None))


def expm1_pred(p: np.ndarray) -> np.ndarray:
    return np.clip(np.expm1(p), 0, None)
