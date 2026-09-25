#!/usr/bin/env python3
"""Parse EPW ``prtgkk`` output into JSON.

Reads the Wannier-interpolated electron-phonon matrix elements |g(k, q, nu)|
printed by EPW (step 8). The table rows share the ph.x prt schema; only the
q/k header format differs. For the diagonal (ibnd == jbnd) band closest to the
Fermi level (the CBM), it reports the rank-1 |g| and the gauge-invariant sum of
|g|^2 over the four highest-omega modes, so the interpolated values can be
compared directly against the step-5 DFPT reference.

Usage:
    python parse_epw_prtgkk.py <epw.out> [--fermi <eV>]

Requirements:
    - Conda environment: base-agent
    - Required packages: none (Python standard library only)
"""

from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import json
import re
from collections import defaultdict

from src.utils.analysis.epw import ROW_RE as ROW_RE


from src.utils.analysis.epw import parse_rows as parse_rows


def summarize(text: str, fermi: float, n_modes: int = 9) -> dict:
    """Summarize the CBM electron-phonon coupling from an EPW prtgkk output.

    Args:
        text: Full contents of the EPW output file.
        fermi: Fermi energy in eV used to select the CBM band.
        n_modes: Number of phonon modes per (k, q) block (3 * natoms).

    Returns:
        Dictionary with the q coordinate, the rank-1 |g| and mode, and the
        gauge-invariant sum of |g|^2 over the four highest-omega modes.
    """
    q = re.search(r"iq\s*=\s*\d+\s+coord\.:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", text)
    rows = parse_rows(text)

    result: dict = {
        "q_cartesian": [float(q.group(i)) for i in (1, 2, 3)] if q else None,
        "fermi_eV_used": fermi,
        "n_rows": len(rows),
    }

    by_ibnd: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        if r["ibnd"] == r["jbnd"]:
            by_ibnd[r["ibnd"]].append(r)
    if by_ibnd:
        cbm = min(by_ibnd, key=lambda b: abs(by_ibnd[b][0]["enk"] - fermi))
        block = sorted(by_ibnd[cbm][:n_modes], key=lambda r: r["imode"])
        gs = [r["g_meV"] for r in block]
        om = [r["omega_meV"] for r in block]
        imax = max(range(len(gs)), key=lambda i: gs[i])
        top4 = sorted(range(len(om)), key=lambda i: om[i], reverse=True)[:4]
        result["cbm_ibnd"] = cbm
        result["cbm_enk_eV"] = block[0]["enk"]
        result["rank1_mode_index"] = block[imax]["imode"]
        result["rank1_omega_meV"] = om[imax]
        result["rank1_g_meV"] = gs[imax]
        result["sum_g2_LO_TO_quartet_meV2"] = sum(gs[i] ** 2 for i in top4)
    return result


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(description="Parse EPW prtgkk output into JSON.")
    parser.add_argument("epw_out", help="Path to the EPW prtgkk output file")
    parser.add_argument(
        "--fermi",
        type=float,
        default=-5.655843,
        help="Fermi energy in eV used to select the CBM band (default: ZrS2 value)",
    )
    parser.add_argument(
        "--n-modes",
        type=int,
        default=9,
        help="Number of phonon modes per (k, q) block, i.e. 3 * natoms (default: 9)",
    )
    args = parser.parse_args()

    with open(args.epw_out) as fh:
        text = fh.read()
    print(json.dumps(summarize(text, args.fermi, args.n_modes), indent=2))


if __name__ == "__main__":
    main()
