"""Existing mat-random-structure-search/scripts/generate_random_structures.py operations shared with installed analysis."""

import warnings
import numpy as np
from pymatgen.core import Composition, Lattice, Structure, Element
from pymatgen.symmetry.groups import SpaceGroup


def estimate_volume_per_atom(elements: list[Element]) -> float:
    """
    Estimate a reasonable volume per atom based on element radii.

    Args:
        elements: List of pymatgen Element objects.

    Returns:
        Estimated volume per atom in Å³.
    """
    radii = []
    for el in elements:
        r = el.atomic_radius
        if r is None:
            r = 1.5  # fallback
        radii.append(float(r))
    avg_radius = np.mean(radii)
    # Approximate: V ≈ (2*r)^3 * packing_factor
    # Use packing factor ~0.7 for random structures
    return (2.0 * avg_radius) ** 3 * 0.7


def get_min_distance(el1: Element, el2: Element) -> float:
    """
    Get minimum allowed distance between two elements.

    Uses sum of covalent radii * 0.7 as a lower bound.

    Args:
        el1, el2: pymatgen Element objects.

    Returns:
        Minimum distance in Å.
    """
    r1 = float(el1.atomic_radius or 1.5)
    r2 = float(el2.atomic_radius or 1.5)
    return (r1 + r2) * 0.7


def check_min_distances(structure: Structure, scale: float = 0.7) -> bool:
    """
    Check that all interatomic distances are above minimum thresholds.

    Args:
        structure: pymatgen Structure to check.
        scale: Scale factor for minimum distance (default: 0.7 * sum of radii).

    Returns:
        True if all distances are acceptable, False otherwise.
    """
    for i, site_i in enumerate(structure):
        for j, site_j in enumerate(structure):
            if j <= i:
                continue
            d = structure.get_distance(i, j)
            min_d = get_min_distance(
                site_i.specie
                if hasattr(site_i.specie, "symbol")
                else Element(str(site_i.specie)),
                site_j.specie
                if hasattr(site_j.specie, "symbol")
                else Element(str(site_j.specie)),
            )
            if d < min_d:
                return False
    return True


def generate_random_structure(
    composition: Composition,
    spacegroup: int,
    volume_scale: float = 1.0,
    max_attempts: int = 50,
) -> Structure | None:
    """
    Generate a random crystal structure for a composition in a given space group.

    Uses pymatgen's Structure.from_spacegroup() with randomized lattice
    parameters and random fractional coordinates for Wyckoff positions.

    Args:
        composition: Target pymatgen Composition.
        spacegroup: Space group number (1-230).
        volume_scale: Scale factor for estimated volume (randomized around this).
        max_attempts: Number of random attempts before giving up.

    Returns:
        A pymatgen Structure, or None if generation failed.
    """
    elements = list(composition.as_dict().keys())
    amounts = [int(composition.as_dict()[el]) for el in elements]
    num_atoms = sum(amounts)

    # Estimate target volume
    element_objs = [Element(el) for el in elements]
    vol_per_atom = estimate_volume_per_atom(element_objs)
    target_volume = vol_per_atom * num_atoms * volume_scale

    sg = SpaceGroup.from_int_number(spacegroup)
    crystal_system = sg.crystal_system

    for _attempt in range(max_attempts):
        # Generate random lattice parameters
        a, b, c, alpha, beta, gamma = _random_lattice_params(
            crystal_system, target_volume
        )

        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Generate random fractional coordinates
        species = []
        coords = []
        for el, amt in zip(elements, amounts):
            for _ in range(amt):
                species.append(el)
                coords.append(np.random.rand(3).tolist())

        # Create structure
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = Structure(
                lattice,
                species,
                coords,
                coords_are_cartesian=False,
            )

        # Check minimum distances
        if check_min_distances(structure):
            return structure

    return None


def _random_lattice_params(
    crystal_system: str,
    target_volume: float,
) -> tuple[float, float, float, float, float, float]:
    """
    Generate random lattice parameters consistent with a crystal system.

    Args:
        crystal_system: One of "cubic", "hexagonal", "trigonal",
                       "tetragonal", "orthorhombic", "monoclinic", "triclinic".
        target_volume: Target unit cell volume in ų.

    Returns:
        Tuple of (a, b, c, alpha, beta, gamma).
    """
    # Random aspect ratios
    r1 = np.random.uniform(0.5, 2.0)
    r2 = np.random.uniform(0.5, 2.0)

    if crystal_system == "cubic":
        a = target_volume ** (1.0 / 3.0)
        return (a, a, a, 90, 90, 90)

    elif crystal_system == "hexagonal":
        # a = b, gamma = 120
        a = (target_volume / (r1 * np.sin(np.radians(60)))) ** (1.0 / 3.0)
        c = a * r1
        return (a, a, c, 90, 90, 120)

    elif crystal_system == "trigonal":
        a = (target_volume / (r1 * np.sin(np.radians(60)))) ** (1.0 / 3.0)
        c = a * r1
        return (a, a, c, 90, 90, 120)

    elif crystal_system == "tetragonal":
        a = (target_volume / r1) ** (1.0 / 3.0)
        c = a * r1
        return (a, a, c, 90, 90, 90)

    elif crystal_system == "orthorhombic":
        a = (target_volume / (r1 * r2)) ** (1.0 / 3.0)
        b = a * r1
        c = a * r2
        return (a, b, c, 90, 90, 90)

    elif crystal_system == "monoclinic":
        beta = np.random.uniform(90, 130)
        a = (target_volume / (r1 * r2 * np.sin(np.radians(beta)))) ** (1.0 / 3.0)
        b = a * r1
        c = a * r2
        return (a, b, c, 90, beta, 90)

    else:  # triclinic
        alpha = np.random.uniform(70, 110)
        beta = np.random.uniform(70, 110)
        gamma = np.random.uniform(70, 110)
        cos_a, cos_b, cos_g = (
            np.cos(np.radians(alpha)),
            np.cos(np.radians(beta)),
            np.cos(np.radians(gamma)),
        )
        sin_g = np.sin(np.radians(gamma))
        vol_factor = np.sqrt(
            1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g
        )
        a = (target_volume / (r1 * r2 * sin_g * vol_factor)) ** (1.0 / 3.0)
        b = a * r1
        c = a * r2
        return (a, b, c, alpha, beta, gamma)
