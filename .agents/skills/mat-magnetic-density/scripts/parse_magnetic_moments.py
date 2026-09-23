"""
Parse magnetic moments from atomate2 VASP calculation results.

This script extracts site-resolved magnetic moments, total magnetization,
and analyzes the magnetic ordering pattern from DFT calculation results.

Usage:
    python parse_magnetic_moments.py results.json --output analysis.json

Requirements:
    - Conda environment: base-agent
    - Required packages: pymatgen, json
"""

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import json
from typing import Dict, Any


from src.utils.analysis.magnetism import parse_magnetic_moments as parse_magnetic_moments


from src.utils.analysis.magnetism import classify_magnetic_ordering as classify_magnetic_ordering


def format_output(analysis: Dict[str, Any]) -> str:
    """
    Format the analysis results as a human-readable string.

    Args:
        analysis: Analysis results dictionary

    Returns:
        Formatted string representation
    """
    lines = []
    lines.append("=" * 60)
    lines.append("MAGNETIC MOMENT ANALYSIS")
    lines.append("=" * 60)

    if analysis["structure_formula"]:
        lines.append(f"\nFormula: {analysis['structure_formula']}")

    if analysis["total_magnetization"] is not None:
        lines.append(f"\nTotal Magnetization: {analysis['total_magnetization']:.3f} μB")

    lines.append(f"Magnetic Ordering: {analysis['magnetic_ordering']}")

    if analysis["site_moments"]:
        lines.append(f"\nNumber of sites: {len(analysis['site_moments'])}")
        lines.append("\nSite-resolved magnetic moments (μB):")
        lines.append("-" * 40)
        lines.append(f"{'Site':<8} {'Element':<10} {'Moment (μB)':<15}")
        lines.append("-" * 40)

        for i, (moment, species) in enumerate(
            zip(
                analysis["site_moments"],
                analysis["species"]
                if analysis["species"]
                else ["?"] * len(analysis["site_moments"]),
            ),
            1,
        ):
            lines.append(f"{i:<8} {species:<10} {moment:>10.3f}")

    lines.append("=" * 60)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Parse magnetic moments from atomate2 VASP calculation results"
    )
    parser.add_argument(
        "input_file", type=str, help="Path to JSON file containing atomate2 results"
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Path to save analysis results as JSON (optional)",
    )

    args = parser.parse_args()

    # Load input data
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)

    with open(input_path, "r") as f:
        results_data = json.load(f)

    # Parse magnetic moments
    analysis = parse_magnetic_moments(results_data)

    # Print formatted output
    print(format_output(analysis))

    # Save to JSON if requested
    if args.output:
        output_path = Path(args.output)
        with open(output_path, "w") as f:
            json.dump(analysis, f, indent=2)
        print(f"\nAnalysis saved to: {output_path}")

    # Save input configs for reproducibility
    from src.utils.config_utils import save_skill_inputs

    save_skill_inputs(args, args.output)


if __name__ == "__main__":
    main()
