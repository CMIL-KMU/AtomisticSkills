"""Parse an EDI ``prefix_inv_tau.dat`` file of state-resolved scattering rates.

EDI's transport pass writes ``prefix_inv_tau.dat`` with one header comment line
and one row per (k-point, band) state that has a nonzero rate, as emitted by
``transport.f90`` (header at line 634, data at lines 637-643):

    ik  ibnd  E(eV)  inv_tau_SERTA(Ry)  inv_tau_MRTA(Ry)  tau_SERTA(fs)  tau_MRTA(fs)

This script loads those rows, reports simple statistics (state count, energy
range, mean SERTA/MRTA lifetimes), and writes a JSON summary plus an
``input_configs.yaml``.

Usage:
    python parse_inv_tau.py mos2_inv_tau.dat --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

import argparse
import json
import os

import yaml

from atomistic_analysis.edi_lifetime import COLUMNS as COLUMNS


from atomistic_analysis.edi_lifetime import parse_inv_tau as parse_inv_tau


from atomistic_analysis.edi_lifetime import summarize as summarize


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse an EDI inv_tau.dat file of state-resolved lifetimes."
    )
    parser.add_argument("inv_tau_dat", help="Path to prefix_inv_tau.dat")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    rows = parse_inv_tau(args.inv_tau_dat)
    stats = summarize(rows)
    os.makedirs(args.output_dir, exist_ok=True)

    for key, value in stats.items():
        print(f"{key:>20}: {value}")

    with open(os.path.join(args.output_dir, "inv_tau_summary.json"), "w") as handle:
        json.dump({"statistics": stats, "states": rows}, handle, indent=2)

    configs = {
        "input_file": os.path.abspath(args.inv_tau_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "columns": COLUMNS,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
