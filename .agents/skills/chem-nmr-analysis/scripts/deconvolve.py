#!/usr/bin/env python3
"""
Deconvolve an NMR mixture against component spectra using Wasserstein (optimal
transport) distance, solved via scipy.optimize.linprog.

Algorithm adapted from Magnetstein/Masserstein (Domzal et al., Anal. Chem. 2024;
Ciach et al., Rapid Commun. Mass Spectrom. 2020).  The dual LP formulation is:

    maximize   v^T z
    subject to T_j^T z <= 0   for each component j
               |z_i - z_{i+1}| <= l_i   (Lipschitz / 1-Wasserstein)
               z_i <= kappa              (denoising penalty)

Component proportions are the dual variables (shadow prices) of the T_j
constraints.  No PuLP, CBC, or external solver required -- uses SciPy's
built-in HiGHS backend.

Usage:
  # Env: nmr-agent
  python deconvolve.py crude.csv ref_a.csv ref_b.csv \
      --protons 18 18 --names borneol isoborneol --baseline-correct --json

Requirements: numpy, scipy (>= 1.7), matplotlib (optional, for --plot)
"""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import json
import pathlib
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Spectrum I/O helpers
# ---------------------------------------------------------------------------


def detect_delim(path: str, default: str = ",") -> str:
    p = pathlib.Path(path)
    if p.suffix.lower() in (".tsv", ".xy"):
        return "\t"
    try:
        with open(path, "r", errors="ignore") as f:
            line = f.readline()
        return "\t" if "\t" in line else default
    except Exception:
        return default


def load_xy(
    path: str, delimiter: Optional[str] = None, mnova: bool = False
) -> np.ndarray:
    if delimiter is None:
        delimiter = "\t" if mnova else detect_delim(path)
    try:
        arr = np.loadtxt(path, delimiter=delimiter, usecols=[0, 1])
    except ValueError:
        arr = np.loadtxt(path, delimiter=delimiter, usecols=[0, 1], skiprows=1)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"{path}: expected two numeric columns (ppm, intensity)")
    return arr


def baseline_correct(arr: np.ndarray) -> np.ndarray:
    corrected = arr.copy()
    corrected[:, 1] -= arr[:, 1].min()
    return corrected


# ---------------------------------------------------------------------------
# Wasserstein LP deconvolution (scipy port of dualdeconv2)
# ---------------------------------------------------------------------------


from atomistic_analysis.nmr import _merge_axes as _merge_axes


from atomistic_analysis.nmr import _intensities_on_axis as _intensities_on_axis


from atomistic_analysis.nmr import _normalize_confs as _normalize_confs


from atomistic_analysis.nmr import wasserstein_deconvolve as wasserstein_deconvolve


from atomistic_analysis.nmr import deconvolve_spectra as deconvolve_spectra


# ---------------------------------------------------------------------------
# Plotting (delegated to plot.py)
# ---------------------------------------------------------------------------


def _save_plot(
    out_path: str,
    mix_arr: np.ndarray,
    comp_arrays: list,
    names: list,
    props: list,
    wd: float,
) -> None:
    from plot import plot_deconvolution

    plot_deconvolution(mix_arr, comp_arrays, names, props, wd, out_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Estimate component proportions in an NMR mixture."
    )
    ap.add_argument(
        "mixture", help="CSV/TSV file with mixture: columns [ppm, intensity]"
    )
    ap.add_argument(
        "components",
        nargs="+",
        help="CSV/TSV files for components: columns [ppm, intensity]",
    )
    ap.add_argument(
        "--protons",
        type=int,
        nargs="+",
        help="Proton counts for each component (e.g. 18 18)",
    )
    ap.add_argument(
        "--names", nargs="+", help="Names for components (e.g. borneol isoborneol)"
    )
    ap.add_argument(
        "--kappa", type=float, default=0.25, help="Denoising penalty (default: 0.25)"
    )
    ap.add_argument(
        "--mnova",
        action="store_true",
        help="Treat inputs as Mnova TSV (delimiter='\\t')",
    )
    ap.add_argument(
        "--baseline-correct",
        action="store_true",
        help="Shift each spectrum so its minimum intensity becomes 0.",
    )
    ap.add_argument(
        "--plot",
        metavar="FILE",
        default=None,
        help="Save a deconvolution plot to FILE (e.g. result.png).",
    )
    ap.add_argument(
        "--json", action="store_true", help="Print JSON result line as well"
    )
    ap.add_argument("--quiet", action="store_true", help="Suppress verbose output")
    args = ap.parse_args()

    # Load data
    mix_arr = load_xy(args.mixture, mnova=args.mnova)
    comp_arrays = [load_xy(p, mnova=args.mnova) for p in args.components]

    if args.baseline_correct:
        if not args.quiet:
            print(
                "Baseline correction: shifting each spectrum so its minimum intensity = 0."
            )
        mix_arr = baseline_correct(mix_arr)
        comp_arrays = [baseline_correct(a) for a in comp_arrays]

    n = len(comp_arrays)
    names = (
        args.names
        if args.names and len(args.names) == n
        else [f"comp{i}" for i in range(n)]
    )
    if args.names and len(args.names) != n:
        print(
            "WARNING: --names length does not match number of components; using default names.",
            file=sys.stderr,
        )

    if args.protons and len(args.protons) != n:
        print(
            "ERROR: --protons length must equal number of components.", file=sys.stderr
        )
        sys.exit(2)
    protons = args.protons if args.protons else [1] * n
    if not args.protons:
        print(
            "NOTE: No --protons provided; assuming 1 for each component.",
            file=sys.stderr,
        )

    # Deconvolve
    result = deconvolve_spectra(mix_arr, comp_arrays, protons, kappa=args.kappa)
    props = result["proportions"]
    wd = result["wasserstein_distance"]

    # Output
    print("\nEstimated proportions:")
    width = max(len(s) for s in names) + 2
    for name, val in zip(names, props):
        print(f"  {name.ljust(width)} {val:.6f}")
    print(f"\nWasserstein distance: {wd:.12f}")

    if args.json:
        out = {"proportions": dict(zip(names, props)), "Wasserstein distance": wd}
        print("\nJSON:", json.dumps(out))

    if args.plot:
        _save_plot(args.plot, mix_arr, comp_arrays, names, props, wd)
        if not args.quiet:
            print(f"\nPlot saved -> {args.plot}")

    # Save input configs for reproducibility
    from src.utils.config_utils import save_skill_inputs

    save_skill_inputs(args, args.plot)


if __name__ == "__main__":
    main()
