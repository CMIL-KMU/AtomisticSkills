"""Reusable functions relocated from chem-dft-orca-advanced-calculation/scripts/parse_orca_output.py without formula changes."""

import re

HARTREE_TO_EV = 27.211386245988


def parse_energy(content: str) -> dict:
    """Extract final single-point energy and component energies."""
    result = {}

    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)", content)
    if matches:
        e_hartree = float(matches[-1])
        result["final_energy_hartree"] = e_hartree
        result["final_energy_eV"] = e_hartree * HARTREE_TO_EV

    nuc_match = re.search(r"Nuclear Repulsion\s+:\s+([-\d.]+)\s+Eh", content)
    if nuc_match:
        result["nuclear_repulsion_hartree"] = float(nuc_match.group(1))

    total_match = re.search(r"Total Energy\s+:\s+([-\d.]+)\s+Eh", content)
    if total_match:
        result["total_energy_hartree"] = float(total_match.group(1))

    disp_match = re.search(r"Dispersion correction\s+([-\d.]+)", content)
    if disp_match:
        result["dispersion_correction_hartree"] = float(disp_match.group(1))

    return result
