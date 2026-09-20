"""Reusable functions relocated from ml-mlip-benchmark/scripts/run_benchmark.py without formula changes."""

import numpy as np


def evaluate_metrics(preds, targets):
    p = np.array(preds).flatten()
    t = np.array(targets).flatten()
    if len(p) == 0 or len(t) == 0:
        return {"mae": 0.0, "rmse": 0.0}
    mae = np.mean(np.abs(p - t))
    rmse = np.sqrt(np.mean((p - t) ** 2))
    return {"mae": float(mae), "rmse": float(rmse)}
