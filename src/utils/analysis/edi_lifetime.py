"""Existing mat-edi-mobility/scripts/parse_inv_tau.py operations shared with installed analysis."""

COLUMNS = [
    "ik",
    "ibnd",
    "E_eV",
    "inv_tau_SERTA_Ry",
    "inv_tau_MRTA_Ry",
    "tau_SERTA_fs",
    "tau_MRTA_fs",
]


def parse_inv_tau(path: str) -> list[dict[str, float]]:
    """Parse an EDI inv_tau.dat file.

    Args:
        path: Path to a ``prefix_inv_tau.dat`` file.

    Returns:
        A list of dicts keyed by COLUMNS. ``ik`` and ``ibnd`` are ints; the
        remaining fields are floats.
    """
    rows: list[dict[str, float]] = []
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < len(COLUMNS):
                continue
            row: dict[str, float] = {}
            for i, col in enumerate(COLUMNS):
                row[col] = int(fields[i]) if col in ("ik", "ibnd") else float(fields[i])
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, float]:
    """Compute summary statistics over parsed inverse-lifetime rows.

    Args:
        rows: Output of :func:`parse_inv_tau`.

    Returns:
        Dictionary with keys: 'n_states', 'E_min_eV', 'E_max_eV',
        'mean_tau_SERTA_fs', 'mean_tau_MRTA_fs'.
    """
    n = len(rows)
    energies = [r["E_eV"] for r in rows]
    return {
        "n_states": n,
        "E_min_eV": min(energies),
        "E_max_eV": max(energies),
        "mean_tau_SERTA_fs": sum(r["tau_SERTA_fs"] for r in rows) / n,
        "mean_tau_MRTA_fs": sum(r["tau_MRTA_fs"] for r in rows) / n,
    }
