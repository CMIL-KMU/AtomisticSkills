"""Existing mat-magnetic-density/scripts/parse_magnetic_moments.py operations shared with installed analysis."""

from typing import Dict, List, Any


def parse_magnetic_moments(results_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract and analyze magnetic moments from atomate2 results.

    Args:
        results_data: Dictionary containing atomate2 calculation results

    Returns:
        Dictionary containing magnetic analysis with keys:
        - 'total_magnetization': Total magnetic moment (μB)
        - 'site_moments': List of magnetic moments per site
        - 'magnetic_ordering': Classification of ordering type
        - 'structure': Structure information
    """
    analysis = {
        "total_magnetization": None,
        "site_moments": [],
        "species": [],
        "magnetic_ordering": "unknown",
        "structure_formula": None,
    }

    # Extract data from results
    if "data" in results_data:
        data = results_data["data"]

        # Handle both list and single structure cases
        if isinstance(data, list) and len(data) > 0:
            calc_data = data[0]
        else:
            calc_data = data

        # Extract structure information
        if "structure" in calc_data:
            structure_data = calc_data["structure"]
            if "sites" in structure_data:
                for site in structure_data["sites"]:
                    if "species" in site:
                        species_list = site["species"]
                        if species_list:
                            element = species_list[0].get("element", "Unknown")
                            analysis["species"].append(element)

            # Try to get formula
            if "composition" in structure_data:
                analysis["structure_formula"] = structure_data.get("formula", "Unknown")

        # Extract magnetic moments - check various possible locations
        if "magmom" in calc_data:
            site_moments = calc_data["magmom"]
            if isinstance(site_moments, (list, tuple)):
                analysis["site_moments"] = [float(m) for m in site_moments]
        elif "output" in calc_data and "magnetic_moments" in calc_data["output"]:
            site_moments = calc_data["output"]["magnetic_moments"]
            if isinstance(site_moments, (list, tuple)):
                analysis["site_moments"] = [float(m) for m in site_moments]

        # Calculate total magnetization if site moments exist
        if analysis["site_moments"]:
            analysis["total_magnetization"] = sum(analysis["site_moments"])

            # Classify magnetic ordering
            analysis["magnetic_ordering"] = classify_magnetic_ordering(
                analysis["site_moments"]
            )

    return analysis


def classify_magnetic_ordering(moments: List[float], threshold: float = 0.1) -> str:
    """
    Classify the type of magnetic ordering based on site moments.

    Args:
        moments: List of magnetic moments per site
        threshold: Threshold below which moments are considered zero (μB)

    Returns:
        String describing the magnetic ordering type
    """
    if not moments:
        return "unknown"

    # Count positive, negative, and near-zero moments
    positive = sum(1 for m in moments if m > threshold)
    negative = sum(1 for m in moments if m < -threshold)
    zero = sum(1 for m in moments if abs(m) <= threshold)

    total_mag = abs(sum(moments))

    if zero == len(moments):
        return "non-magnetic"
    elif (negative == 0 and positive > 0) or (positive == 0 and negative > 0):
        return "ferromagnetic"
    elif positive > 0 and negative > 0:
        if total_mag < threshold:
            return "antiferromagnetic"
        else:
            return "ferrimagnetic"
    else:
        return "unknown"
