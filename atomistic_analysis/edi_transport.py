"""Existing mat-edi-mobility/scripts/parse_transport.py operations shared with installed analysis."""

COLUMNS = ["T_K", "mu_SERTA_xx", "mu_MRTA_xx", "mu_SERTA_yy", "mu_MRTA_yy"]


def parse_transport(path: str) -> tuple[list[dict[str, float]], dict[str, str]]:
    """Parse an EDI transport.dat file.

    Args:
        path: Path to a ``prefix_transport.dat`` file.

    Returns:
        A tuple ``(rows, meta)`` where ``rows`` is a list of dicts keyed by
        COLUMNS (all cm^2/Vs except T_K in kelvin) and ``meta`` holds the
        comment-header strings (grid, window) for provenance.
    """
    rows: list[dict[str, float]] = []
    meta: dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if "Grid" in stripped:
                    meta["grid"] = stripped.lstrip("# ").strip()
                elif "Window" in stripped:
                    meta["window"] = stripped.lstrip("# ").strip()
                continue
            fields = stripped.split()
            if len(fields) < len(COLUMNS):
                continue
            rows.append({col: float(fields[i]) for i, col in enumerate(COLUMNS)})
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows, meta
