"""Reusable functions relocated from chem-spectrum-matcher/scripts/match_spectrum.py without formula changes."""

import numpy as np


def normalize(y: np.ndarray) -> np.ndarray:
    """Min-max normalize to [0, 1]. Returns zeros if flat."""
    y = y - y.min()
    m = y.max()
    return y / m if m > 0 else y


def similarity_l2(y1: np.ndarray, y2: np.ndarray) -> float:
    """1 - normalized L2 distance. Range [0, 1], higher = more similar."""
    dist = np.sqrt(np.mean((y1 - y2) ** 2))
    return float(max(0.0, 1.0 - dist))


def similarity_cosine(y1: np.ndarray, y2: np.ndarray) -> float:
    """Cosine similarity. Range [0, 1]."""
    n1, n2 = np.linalg.norm(y1), np.linalg.norm(y2)
    if n1 == 0 or n2 == 0:
        return 0.0
    return float(np.clip(np.dot(y1, y2) / (n1 * n2), 0.0, 1.0))


def similarity_wasserstein(y1: np.ndarray, y2: np.ndarray) -> float:
    """
    1 - normalized Wasserstein-1 distance between two normalized distributions.

    Treats spectra as probability distributions. Range [0, 1].
    """
    from scipy.stats import wasserstein_distance

    # Normalize to probability distributions
    s1, s2 = y1.sum(), y2.sum()
    if s1 == 0 or s2 == 0:
        return 0.0
    dist = wasserstein_distance(
        np.arange(len(y1)), np.arange(len(y2)), y1 / s1, y2 / s2
    )
    # Normalize by max possible distance (full width)
    max_dist = len(y1)
    return float(max(0.0, 1.0 - dist / max_dist))
