"""Existing scalar skill functions, relocated without changing their formulae."""

import json
from pathlib import Path
from typing import Dict, Any

EV_PER_A2_TO_J_PER_M2 = 16.0218


def calculate_voltage(
    e_full: float,
    e_empty: float,
    e_metal: float,
    n_metal: int,
    n_ions: int,
    metal_symbol: str = None,
    output_file: str = None,
) -> Dict[str, Any]:
    """
    Calculate the average intercalation voltage.

    The voltage is calculated using:
    V = -(E_full - E_empty - n * μ_metal) / n

    where μ_metal = E_metal / n_metal_atoms

    Args:
        e_full: Total energy of fully intercalated structure (eV)
        e_empty: Total energy of de-intercalated structure (eV)
        e_metal: Total energy of bulk metal structure (eV)
        n_metal: Number of metal atoms in the bulk metal structure
        n_ions: Number of intercalated ions (n_full - n_empty)
        metal_symbol: Optional symbol of intercalating ion for documentation
        output_file: Optional path to save results as JSON

    Returns:
        Dictionary with voltage and energy data
    """
    # Calculate chemical potential of metal
    mu_metal = e_metal / n_metal

    # Calculate voltage
    # V = -(E_full - E_empty - n*μ_metal) / n
    voltage = -(e_full - e_empty - n_ions * mu_metal) / n_ions

    results = {
        "voltage_V": voltage,
        "n_ions": n_ions,
        "E_full_eV": e_full,
        "E_empty_eV": e_empty,
        "E_metal_total_eV": e_metal,
        "n_metal_atoms": n_metal,
        "mu_metal_eV": mu_metal,
        "energy_difference_eV": e_full - e_empty,
        "metal_contribution_eV": n_ions * mu_metal,
    }

    if metal_symbol:
        results["metal_symbol"] = metal_symbol

    # Print results
    print("\n" + "=" * 60)
    print("INTERCALATION VOLTAGE CALCULATION")
    print("=" * 60)
    print(f"Average Voltage: {voltage:.4f} V")
    print("\nEnergies:")
    print(f"  E(full):        {e_full:.6f} eV")
    print(f"  E(empty):       {e_empty:.6f} eV")
    print(f"  E(metal):       {e_metal:.6f} eV ({n_metal} atoms)")
    print(f"  μ(metal):       {mu_metal:.6f} eV/atom")
    print("\nIntercalation:")
    print(f"  Number of ions: {n_ions}")
    print(f"  ΔE:             {e_full - e_empty:.6f} eV")
    print(f"  n×μ(metal):     {n_ions * mu_metal:.6f} eV")
    print("=" * 60 + "\n")

    # Save results if output file specified
    if output_file:
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to: {output_path}\n")

    return results


def compute_gb_energy(
    e_gb: float,
    n_atoms: int,
    e_bulk_per_atom: float,
    area_A2: float,
) -> float:
    """Return grain boundary energy in J/m²."""
    delta_e = e_gb - n_atoms * e_bulk_per_atom  # eV
    gamma = delta_e / (2.0 * area_A2)  # eV/Å²
    return gamma * EV_PER_A2_TO_J_PER_M2  # J/m²
