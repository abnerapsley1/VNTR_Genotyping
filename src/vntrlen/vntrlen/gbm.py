"""Gradient-boosted trees (XGBoost), saved as portable JSON boosters.

Training settings match the shipped models: leaf-wise growth, 63 leaves, L2 1.0,
255 bins, learning rate 0.05, up to 600 rounds with early stopping on a random 10%.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import xgboost as xgb


def _device() -> str:
    try:
        xgb.train({"device": "cuda", "tree_method": "hist"}, xgb.DMatrix(np.zeros((2, 1)), label=[0, 1]), 1)
        return "cuda"
    except xgb.core.XGBoostError:
        return "cpu"


DEVICE = None


def device() -> str:
    global DEVICE
    if DEVICE is None:
        DEVICE = _device()
    return DEVICE


def train(X: np.ndarray, y: np.ndarray, quantile: float | None = None, seed: int = 20260831,
          n_estimators: int = 600, learning_rate: float = 0.05) -> xgb.Booster:
    params = dict(objective="reg:quantileerror" if quantile is not None else "reg:squarederror",
                  learning_rate=learning_rate, grow_policy="lossguide", max_leaves=63, max_depth=0,
                  reg_lambda=1.0, max_bin=255, tree_method="hist", device=device(), seed=seed)
    if quantile is not None:
        params["quantile_alpha"] = quantile
    idx = np.random.default_rng(seed).permutation(len(y))
    nv = max(1, int(0.1 * len(y)))
    va, tr = idx[:nv], idx[nv:]
    dtr = xgb.QuantileDMatrix(X[tr], y[tr], max_bin=255)
    dva = xgb.DMatrix(X[va], y[va])
    return xgb.train(params, dtr, n_estimators, evals=[(dva, "val")], early_stopping_rounds=30,
                     verbose_eval=False)


def predict(b: xgb.Booster, X: np.ndarray) -> np.ndarray:
    b.set_param({"device": device()})
    rng = (0, b.best_iteration + 1) if b.attributes().get("best_iteration") else (0, 0)
    with warnings.catch_warnings():       # host-memory input on a GPU booster: correct, just slower
        warnings.simplefilter("ignore", UserWarning)
        return b.inplace_predict(np.ascontiguousarray(X, dtype=np.float32), iteration_range=rng)


def save(b: xgb.Booster, path: Path) -> None:
    b.save_model(str(path))


def load(path: Path) -> xgb.Booster:
    b = xgb.Booster()
    b.load_model(str(path))
    return b
