"""
Calculate the average intercalation voltage from relaxed structure energies.

This script takes the energies of the fully intercalated state, de-intercalated state,
and bulk metal to compute the average voltage.

Usage:
    python calculate_voltage.py --e_full -123.45 --e_empty -98.76 --e_metal -1.23 --n_metal 16 --n_ions 4

Requirements:
    - Conda environment: base-agent
    - Required packages: argparse, json, pathlib
"""

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from atomistic_analysis.scalars import calculate_voltage


def main():
    parser = argparse.ArgumentParser(
        description="Calculate average intercalation voltage from energies"
    )
    parser.add_argument(
        "--e_full",
        type=float,
        required=True,
        help="Total energy of fully intercalated structure (eV)",
    )
    parser.add_argument(
        "--e_empty",
        type=float,
        required=True,
        help="Total energy of de-intercalated structure (eV)",
    )
    parser.add_argument(
        "--e_metal",
        type=float,
        required=True,
        help="Total energy of bulk metal structure (eV)",
    )
    parser.add_argument(
        "--n_metal",
        type=int,
        required=True,
        help="Number of metal atoms in the bulk metal structure",
    )
    parser.add_argument(
        "--n_ions",
        type=int,
        required=True,
        help="Number of intercalated ions (difference between full and empty)",
    )
    parser.add_argument(
        "--metal",
        type=str,
        default=None,
        help="Symbol of intercalating ion (for documentation)",
    )
    parser.add_argument(
        "--output", type=str, default=None, help="Path to save results as JSON"
    )

    args = parser.parse_args()

    calculate_voltage(
        e_full=args.e_full,
        e_empty=args.e_empty,
        e_metal=args.e_metal,
        n_metal=args.n_metal,
        n_ions=args.n_ions,
        metal_symbol=args.metal,
        output_file=args.output,
    )

    # Save input configs for reproducibility
    from src.utils.config_utils import save_skill_inputs

    save_skill_inputs(args, args.output)


if __name__ == "__main__":
    main()
