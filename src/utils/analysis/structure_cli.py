"""Bounded PNG rendering entry point; accepts explicit coordinates, never file paths."""

import json
import sys


def main():
    from pymatgen.core import Structure
    from src.utils.analysis.structure_viz import structure_3d_custom

    data = json.load(sys.stdin)
    if not 0 < len(data["species"]) <= 500:
        raise ValueError("Preview supports 1–500 sites")
    structure = Structure(data["lattice"], data["species"], data["coords"])
    # No inferred bonds or symmetry standardization: display the recorded coordinates.
    figure = structure_3d_custom(
        structure,
        show_bonds=False,
        atom_size=22,
        site_labels=False,
        show_site_vectors=(),
    )
    figure.update_layout(paper_bgcolor="white", font_color="#001219")
    sys.stdout.buffer.write(
        figure.to_image(format="png", width=1000, height=1000, scale=1)
    )


if __name__ == "__main__":
    main()
