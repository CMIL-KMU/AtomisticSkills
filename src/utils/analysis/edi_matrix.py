"""Existing mat-edi-mobility/scripts/parse_edmat.py operations shared with installed analysis."""


def parse_edmat(path: str) -> tuple[list[dict[str, float]], int]:
    """Parse an EDI matrix-element file, auto-detecting its layout.

    Args:
        path: Path to an EDI ``*_edmat*.dat`` file.

    Returns:
        A tuple ``(records, ncol)`` where ``records`` is a list of normalised
        dicts and ``ncol`` is the detected field count (9, 13, or 15).

    Raises:
        ValueError: If the file has no data rows or an unrecognised layout.
    """
    records: list[dict[str, float]] = []
    ncol = 0
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if ncol == 0:
                ncol = len(fields)
                if ncol not in (9, 13, 15):
                    raise ValueError(
                        f"Unrecognised EDI matrix-element layout: {ncol} columns in {path}"
                    )
            records.append(_parse_row(fields, ncol))
    if not records:
        raise ValueError(f"No data rows found in {path}")
    return records, ncol


def _parse_row(fields: list[str], ncol: int) -> dict[str, float]:
    """Normalise one data row to a common record schema.

    Args:
        fields: Whitespace-split tokens of a single data line.
        ncol: The layout width detected for the file (9, 13, or 15).

    Returns:
        Dict with keys ikf, kf, ibnd, jbnd, m2, reM, imM (+ m2_loc, m2_nl
        for the 15-column direct layout).
    """
    if ncol == 9:
        # ik  kx ky kz  ibnd jbnd  |M|^2  Re  Im
        return {
            "ikf": int(fields[0]),
            "kf": [float(fields[1]), float(fields[2]), float(fields[3])],
            "ibnd": int(fields[4]),
            "jbnd": int(fields[5]),
            "m2": float(fields[6]),
            "reM": float(fields[7]),
            "imM": float(fields[8]),
        }
    # 13- and 15-column layouts share the first 13 fields.
    record = {
        "ikf": int(fields[1]),
        "kf": [float(fields[5]), float(fields[6]), float(fields[7])],
        "ibnd": int(fields[8]),
        "jbnd": int(fields[9]),
        "m2": float(fields[10]),
        "reM": float(fields[11]),
        "imM": float(fields[12]),
    }
    if ncol == 15:
        record["m2_loc"] = float(fields[13])
        record["m2_nl"] = float(fields[14])
    return record


def select_band_pair(
    records: list[dict[str, float]], ibnd: int, jbnd: int
) -> list[dict[str, float]]:
    """Filter records to a single (ibnd, jbnd) band pair, ordered by ikf.

    Args:
        records: Output of :func:`parse_edmat`.
        ibnd: Initial-state band index to keep.
        jbnd: Final-state band index to keep.

    Returns:
        The matching records sorted by final-k index ``ikf``.
    """
    picked = [r for r in records if r["ibnd"] == ibnd and r["jbnd"] == jbnd]
    return sorted(picked, key=lambda r: r["ikf"])
