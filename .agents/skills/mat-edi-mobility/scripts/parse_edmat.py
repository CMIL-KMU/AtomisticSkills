"""Parse EDI electron-defect matrix-element output files ``|M(k_i, k_f)|^2``.

EDI writes several matrix-element text files that all share a fixed-width layout
with ``|M|^2  Re(M)  Im(M)`` as the last physical columns. This parser
auto-detects which of the three EDI layouts a file uses, from the number of
whitespace-separated fields in the first data row:

* 9 fields  -> single-k list (``prefix_edmatf.dat`` diagonal interpolation, or
  ``prefix_edmat_scatter.dat`` from-K interpolation):
  ``ik  kx ky kz  ibnd jbnd  |M|^2  Re(M)  Im(M)``
  (``ed_coarse.f90`` headers at lines 71 and 198).
* 13 fields -> Wannier-interpolated pair list (``prefix_edmat_interp.dat``,
  produced by ``edmat_interp_from_file``):
  ``iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)``
  (``ed_coarse.f90`` header at line 2526, format at line 2551).
* 15 fields -> direct / coarse Bloch pair list (``prefix_edmat_direct.dat`` from
  ``edmat_direct_from_file``, or ``prefix_edmat_bloch.dat``): the 13-field
  layout plus ``|M_loc|^2  |M_nl|^2``
  (``ed_coarse.f90`` headers at lines 2145 and 739, format at line 2405).

Each parsed record is normalised to keys: ``ikf`` (final-k index), ``kf`` (3
floats), ``ibnd``, ``jbnd``, ``m2`` (``|M|^2`` in Ry^2), ``reM``, ``imM``, and,
when present, ``m2_loc`` and ``m2_nl``.

Usage:
    python parse_edmat.py mos2_edmat_interp.dat --ibnd 2 --jbnd 2 --output-dir results/

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


from src.utils.analysis.edi_matrix import parse_edmat as parse_edmat


from src.utils.analysis.edi_matrix import _parse_row as _parse_row


from src.utils.analysis.edi_matrix import select_band_pair as select_band_pair


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse an EDI electron-defect matrix-element file."
    )
    parser.add_argument("edmat_dat", help="Path to an EDI *_edmat*.dat file")
    parser.add_argument("--ibnd", type=int, default=None, help="Keep only this ibnd")
    parser.add_argument("--jbnd", type=int, default=None, help="Keep only this jbnd")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    records, ncol = parse_edmat(args.edmat_dat)
    if args.ibnd is not None and args.jbnd is not None:
        records = select_band_pair(records, args.ibnd, args.jbnd)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Detected {ncol}-column layout; {len(records)} records retained.")
    if records:
        m2_values = [r["m2"] for r in records]
        print(f"max |M|^2 = {max(m2_values):.6e} Ry^2")

    with open(os.path.join(args.output_dir, "edmat_records.json"), "w") as handle:
        json.dump({"ncol": ncol, "records": records}, handle, indent=2)

    configs = {
        "input_file": os.path.abspath(args.edmat_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "ibnd": args.ibnd,
        "jbnd": args.jbnd,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
