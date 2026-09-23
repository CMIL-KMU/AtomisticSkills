"""Existing mat-epw-mobility/scripts/gen_kpoints.py operations shared with installed analysis."""


def uniform_kpoints(n1: int, n2: int, n3: int) -> list[tuple[float, float, float]]:
    """Build the fractional coordinates of a uniform Gamma-centred mesh.

    Args:
        n1: Number of divisions along b1.
        n2: Number of divisions along b2.
        n3: Number of divisions along b3 (use 1 for a 2D slab).

    Returns:
        List of (k1, k2, k3) fractional coordinates.
    """
    return [
        (i / n1, j / n2, k / n3)
        for i in range(n1)
        for j in range(n2)
        for k in range(n3)
    ]


def format_card(points: list[tuple[float, float, float]]) -> str:
    """Format a k-point list as a QE ``K_POINTS crystal`` card.

    Args:
        points: Fractional coordinates from :func:`uniform_kpoints`.

    Returns:
        The card text, including the ``K_POINTS crystal`` header and count line.
    """
    lines = ["K_POINTS crystal", str(len(points))]
    lines += [f"{a:.10f}  {b:.10f}  {c:.10f}  1.0" for a, b, c in points]
    return "\n".join(lines)
