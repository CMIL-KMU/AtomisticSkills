"""Explicit-input analysis CLI. Use the base-agent environment; no project package installation is required.

Use --help for required timing, species, charge and output options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from src.utils.analysis.cli import diffusion_main

if __name__ == "__main__":
    diffusion_main()
