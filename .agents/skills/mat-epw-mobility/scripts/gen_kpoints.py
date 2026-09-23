#!/usr/bin/env python3
"""Generate an explicit uniform k-point list for a QE NSCF ``K_POINTS crystal`` card.

EPW reads back the full uniform mesh, so the NSCF must use an explicit k-list
(``K_POINTS crystal``), never ``automatic`` (step 2). This prints the header
line and the ``n1 * n2 * n3`` fractional points with unit weights, ready to
paste into ``nscf.in``.

Usage:
    python gen_kpoints.py 12 12 1

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


from src.utils.analysis.kpoints import uniform_kpoints as uniform_kpoints


from src.utils.analysis.kpoints import format_card as format_card


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Generate an explicit uniform K_POINTS crystal card for QE NSCF."
    )
    parser.add_argument("n1", type=int, help="Divisions along b1")
    parser.add_argument("n2", type=int, help="Divisions along b2")
    parser.add_argument("n3", type=int, help="Divisions along b3 (1 for a 2D slab)")
    args = parser.parse_args()

    print(format_card(uniform_kpoints(args.n1, args.n2, args.n3)))


if __name__ == "__main__":
    main()
