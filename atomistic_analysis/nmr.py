"""Existing chem-nmr-analysis/scripts/deconvolve.py operations shared with installed analysis."""

import numpy as np
from scipy.optimize import linprog


def _merge_axes(*spectra_confs: list[tuple[float, float]]) -> np.ndarray:
    """Build a sorted, deduplicated common ppm axis from multiple spectra."""
    all_ppm = set()
    for confs in spectra_confs:
        for ppm, _ in confs:
            all_ppm.add(round(ppm, 6))
    return np.array(sorted(all_ppm))


def _intensities_on_axis(
    confs: list[tuple[float, float]], axis: np.ndarray
) -> np.ndarray:
    """Project a (ppm, intensity) list onto a common axis (zero where absent)."""
    lookup = {}
    for ppm, inten in confs:
        key = round(ppm, 6)
        lookup[key] = lookup.get(key, 0.0) + inten
    return np.array([lookup.get(round(a, 6), 0.0) for a in axis])


def _normalize_confs(confs: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Clamp negatives to 0, then normalize intensities to sum to 1."""
    clamped = [(p, max(i, 0.0)) for p, i in confs]
    total = sum(i for _, i in clamped)
    if total <= 0:
        return clamped
    return [(p, i / total) for p, i in clamped]


def wasserstein_deconvolve(
    mixture_confs: list[tuple[float, float]],
    component_confs_list: list[list[tuple[float, float]]],
    kappa: float = 0.25,
) -> dict:
    """
    Estimate component proportions via Wasserstein-distance LP deconvolution.

    Parameters
    ----------
    mixture_confs : list of (ppm, intensity)
        Mixture spectrum, will be normalized internally.
    component_confs_list : list of list of (ppm, intensity)
        Reference spectra, each normalized internally.
    kappa : float
        Denoising penalty (max transport cost).  Default 0.25.

    Returns
    -------
    dict with keys:
        "proportions" : list[float]  -- estimated proportion per component
        "wasserstein_distance" : float  -- objective value (fit quality)
        "noise" : float  -- fraction of signal attributed to noise
    """
    # Normalize
    mix_n = _normalize_confs(mixture_confs)
    comp_n = [_normalize_confs(c) for c in component_confs_list]

    # Build common axis
    axis = _merge_axes(mix_n, *comp_n)
    n = len(axis)
    if n < 2:
        return {
            "proportions": [0.0] * len(comp_n),
            "wasserstein_distance": float("nan"),
            "noise": 1.0,
        }

    # Intensity vectors on common axis
    v = _intensities_on_axis(mix_n, axis)  # mixture
    T = [_intensities_on_axis(c, axis) for c in comp_n]  # components
    k = len(T)

    # Interval lengths between consecutive axis points
    intervals = np.diff(axis)

    # -----------------------------------------------------------------------
    # LP formulation  (scipy minimizes, so we minimize -v^T z)
    #
    # Variables: z_0 .. z_{n-1}
    #
    # Inequality constraints (A_ub @ z <= b_ub):
    #   1. Component constraints:  T_j^T z <= 0    (k constraints)
    #   2. Lipschitz forward:   z_i - z_{i+1} <= l_i   (n-1 constraints)
    #   3. Lipschitz backward:  z_{i+1} - z_i <= l_i   (n-1 constraints)
    #
    # Variable bounds:  z_i <= kappa  (no lower bound)
    # -----------------------------------------------------------------------

    c_obj = -v  # minimize -v^T z  <=>  maximize v^T z

    n_ineq = k + 2 * (n - 1)
    A_ub = np.zeros((n_ineq, n))
    b_ub = np.zeros(n_ineq)

    row = 0
    # Component constraints: T_j^T z <= 0
    for j in range(k):
        A_ub[row, :] = T[j]
        b_ub[row] = 0.0
        row += 1

    # Lipschitz forward: z_i - z_{i+1} <= l_i
    for i in range(n - 1):
        A_ub[row, i] = 1.0
        A_ub[row, i + 1] = -1.0
        b_ub[row] = intervals[i]
        row += 1

    # Lipschitz backward: z_{i+1} - z_i <= l_i
    for i in range(n - 1):
        A_ub[row, i] = -1.0
        A_ub[row, i + 1] = 1.0
        b_ub[row] = intervals[i]
        row += 1

    bounds = [(None, kappa)] * n

    result = linprog(
        c_obj,
        A_ub=A_ub,
        b_ub=b_ub,
        bounds=bounds,
        method="highs",
        options={"presolve": True, "disp": False},
    )

    if not result.success:
        return {
            "proportions": [0.0] * k,
            "wasserstein_distance": float("nan"),
            "noise": 1.0,
        }

    # Extract proportions from dual variables of component constraints
    # result.ineqlin.marginals[j] for the first k constraints
    duals = result.ineqlin.marginals
    proportions = [
        -float(duals[j]) for j in range(k)
    ]  # negated because dual of <= constraint

    # Noise: from reduced costs (dual of upper bound constraints)
    noise = 1.0 - sum(proportions)

    # Wasserstein distance = optimal objective value (negated back)
    wd = -float(result.fun)

    return {
        "proportions": proportions,
        "wasserstein_distance": wd,
        "noise": max(noise, 0.0),
    }


def deconvolve_spectra(
    mix_arr: np.ndarray,
    comp_arrays: list[np.ndarray],
    protons: list[int],
    kappa: float = 0.25,
) -> dict:
    """
    High-level deconvolution: load arrays, run LP, apply proton correction.

    Parameters
    ----------
    mix_arr : ndarray of shape (M, 2)
        Mixture spectrum (ppm, intensity).
    comp_arrays : list of ndarray of shape (N_i, 2)
        Component reference spectra.
    protons : list of int
        Number of 1H protons per component molecule.
    kappa : float
        Denoising penalty.

    Returns
    -------
    dict with "proportions", "wasserstein_distance", "noise"
    """
    mix_confs = list(zip(mix_arr[:, 0].tolist(), mix_arr[:, 1].tolist()))
    comp_confs = [list(zip(a[:, 0].tolist(), a[:, 1].tolist())) for a in comp_arrays]

    raw = wasserstein_deconvolve(mix_confs, comp_confs, kappa=kappa)
    raw_props = raw["proportions"]

    # Proton correction: convert area-proportional to concentration-proportional
    if protons and all(p > 0 for p in protons):
        corrected = [prop / p for prop, p in zip(raw_props, protons)]
        total = sum(corrected)
        if total > 0:
            corrected = [c / total for c in corrected]
        raw["proportions"] = corrected

    return raw
