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
        import xgboost as xgb
        obj = "reg:quantileerror" if quantile is not None else "reg:squarederror"
        kw = dict(objective=obj, n_estimators=n_estimators, learning_rate=learning_rate,
                  max_depth=max_depth or 6, subsample=0.8, colsample_bytree=0.8,
                  random_state=seed, tree_method="hist")
        try:
            kw["device"] = "cuda"
        except Exception:  # noqa: BLE001
            pass
        if quantile is not None:
            kw["quantile_alpha"] = quantile
        return xgb.XGBRegressor(**kw)
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


def fit_predict(Xtr: np.ndarray, ytr: np.ndarray, Xall: np.ndarray, **kw) -> np.ndarray:
    m = make_regressor(**kw)
    m.fit(Xtr, ytr)
    return m.predict(Xall)


def log1p_target(y: np.ndarray) -> np.ndarray:
    return np.log1p(np.clip(y, 0, None))


def expm1_pred(p: np.ndarray) -> np.ndarray:
    return np.clip(np.expm1(p), 0, None)
