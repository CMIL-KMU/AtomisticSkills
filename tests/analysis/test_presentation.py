"""Verify original rendering and upstream style ownership without changing science."""

import matplotlib.pyplot as plt
from matplotlib import rc_params_from_file
from atomistic_analysis.plotting import browser_style, plot_style, style_path


def test_lovelyplots_style_is_shared_and_scoped():
    values = rc_params_from_file(style_path(), use_default_template=False)
    theme = browser_style()
    before = plt.rcParams.copy()
    assert theme["colors"] == values["axes.prop_cycle"].by_key()["color"]
    assert theme["grid"] is False
    with plot_style():
        assert plt.rcParams["axes.prop_cycle"] == values["axes.prop_cycle"]
        assert plt.rcParams["svg.fonttype"] == "none"
    assert plt.rcParams == before


def test_structure_renderer_preserves_recorded_cell():
    from pymatgen.core import Lattice, Structure
    from atomistic_analysis.structure_viz import structure_3d_custom
    from src.utils.structure_viz import structure_3d_custom as legacy

    assert legacy is structure_3d_custom
    structure = Structure(
        Lattice.from_parameters(4, 5, 6, 80, 85, 75),
        ["Li", "Cl"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    original = structure.as_dict()
    figure = structure_3d_custom(structure, show_bonds=False, show_site_vectors=())
    assert structure.as_dict() == original
    assert "Cartesian x, y, z" in figure.layout.title.text
    assert all(
        getattr(figure.layout, s).aspectmode == "data"
        for s in ("scene", "scene2", "scene3", "scene4")
    )
