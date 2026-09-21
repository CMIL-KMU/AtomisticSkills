"""Anvil policy is explicit and does not remove upstream scientific operations."""

import copy
import pytest
from atomistic_analysis import catalog, tools


def test_explicit_selection_exclusions_and_contracts(monkeypatch):
    assert len(catalog.OPERATIONS) == 27
    for name in catalog.CATALOG["excluded_operations"]:
        assert name in tools.SCHEMAS
        with pytest.raises(ValueError):
            catalog.selected(name)
    monkeypatch.setitem(tools.SCHEMAS, "future-unreviewed", {})
    with pytest.raises(ValueError):
        catalog.selected("future-unreviewed")
    contract = catalog.contract("xrd")
    schema = copy.deepcopy(tools.SCHEMAS["xrd"])
    schema["properties"]["wavelength_angstrom"]["maximum"] = 100
    monkeypatch.setitem(tools.SCHEMAS, "xrd", schema)
    assert contract != catalog.contract("xrd")


def test_density_axes_do_not_assume_equal_fractional_spacing():
    axes = catalog.density_axes(
        dict(
            frames=[dict(lattice=[[6, 0, 0], [0, 6, 0], [0, 0, 6]])],
            interval_angstrom=2.5,
        ),
        [3, 3, 3],
    )
    assert axes[0]["labels"] == ["0.0", str(2.5 / 6), str(5 / 6)]
