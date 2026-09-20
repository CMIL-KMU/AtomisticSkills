"""Explicit-input analysis CLI. Install atomistic-analysis[transport] in base-agent.

Use --help for required timing, species, charge and output options.
"""
from atomistic_analysis.cli import diffusion_main

if __name__ == "__main__":
    diffusion_main()
