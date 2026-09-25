"""Existing coordinate box operation shared by recipes and the skill CLI."""

import numpy as np
from typing import Dict


def compute_box(
    positions: np.ndarray,
    padding: float,
    min_size: float,
) -> Dict[str, float]:
    """Compute bounding box center and size from an (N, 3) coordinate array."""
    if len(positions) == 0:
        raise ValueError("No coordinates provided for box computation")

    pos_min = positions.min(axis=0)
    pos_max = positions.max(axis=0)
    center = 0.5 * (pos_min + pos_max)
    size = np.maximum(pos_max - pos_min + 2.0 * padding, min_size)

    return {
        "center_x": round(float(center[0]), 4),
        "center_y": round(float(center[1]), 4),
        "center_z": round(float(center[2]), 4),
        "size_x": round(float(size[0]), 4),
        "size_y": round(float(size[1]), 4),
        "size_z": round(float(size[2]), 4),
        "padding": padding,
        "min_size": min_size,
    }
