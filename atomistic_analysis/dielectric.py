"""Existing mat-dielectric-response/scripts/plot_dielectric.py operations shared with installed analysis."""

from pathlib import Path
import re
import gzip
import numpy as np


def read_text_maybe_gz(path: Path) -> str:
    """
    Read plain-text content from a regular or gzipped file.

    Args:
        path: File path to read.

    Returns:
        Decoded file content.
    """
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="ignore") as handle:
            return handle.read()
    return path.read_text(errors="ignore")


def parse_dielectric_section(
    outcar_path: Path, header: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse a dielectric data block from OUTCAR.

    The relevant tutorial blocks contain rows with energy, Re(epsilon), and
    Im(epsilon). Additional trailing text is ignored.

    Args:
        outcar_path: Path to OUTCAR from an `ALGO = CHI` calculation.
        header: Header line identifying the section to extract.

    Returns:
        Tuple of energies, real part, and imaginary part.
    """
    text = read_text_maybe_gz(outcar_path)
    lines = text.splitlines()

    start_idx = None
    for idx, line in enumerate(lines):
        if header in line:
            start_idx = idx + 1
            break

    if start_idx is None:
        raise ValueError(f"Could not find section '{header}' in {outcar_path}")

    number_pattern = re.compile(
        r"^\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)"
    )

    energies: list[float] = []
    real: list[float] = []
    imag: list[float] = []
    saw_data = False

    for line in lines[start_idx:]:
        match = number_pattern.match(line)
        if match:
            energies.append(float(match.group(1)))
            real.append(float(match.group(2)))
            imag.append(float(match.group(3)))
            saw_data = True
            continue

        if saw_data and line.strip():
            break

    if not energies:
        raise ValueError(
            f"No dielectric data found after section '{header}' in {outcar_path}"
        )

    return np.array(energies), np.array(real), np.array(imag)
