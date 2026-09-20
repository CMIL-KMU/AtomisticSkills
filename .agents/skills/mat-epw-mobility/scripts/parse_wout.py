#!/usr/bin/env python3
"""Parse a Wannier90 ``.wout`` Final State block into JSON.

Extracts the maximally-localized Wannier function (MLWF) spread decomposition
and centres from the last ``Final State`` block. This is the primary quality
gate for the Wannier stage (step 7): a total spread that is large relative to
the gauge-invariant lower bound Omega_I signals image folding in a long-vacuum
2D cell.

Usage:
    python parse_wout.py <prefix.wout>

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


from atomistic_analysis.wannier import parse_wout as parse_wout


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Parse a Wannier90 .wout Final State block into JSON."
    )
    parser.add_argument("wout", help="Path to the Wannier90 .wout file")
    args = parser.parse_args()

    with open(args.wout) as fh:
        text = fh.read()
    print(json.dumps(parse_wout(text), indent=2))


if __name__ == "__main__":
    main()
