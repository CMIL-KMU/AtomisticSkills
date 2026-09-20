"""Existing mat-epw-mobility/scripts/parse_epw_prtgkk.py operations shared with installed analysis."""

import re

ROW_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*$",
    re.MULTILINE,
)


def parse_rows(text: str) -> list[dict]:
    """Extract the |g| table rows from an EPW prtgkk output.

    Args:
        text: Full contents of the EPW output file.

    Returns:
        List of row dictionaries with keys ibnd, jbnd, imode, enk, enkq,
        omega_meV, g_meV.
    """
    rows: list[dict] = []
    for m in ROW_RE.finditer(text):
        rows.append(
            {
                "ibnd": int(m.group(1)),
                "jbnd": int(m.group(2)),
                "imode": int(m.group(3)),
                "enk": float(m.group(4)),
                "enkq": float(m.group(5)),
                "omega_meV": float(m.group(6)),
                "g_meV": float(m.group(7)),
            }
        )
    return rows
