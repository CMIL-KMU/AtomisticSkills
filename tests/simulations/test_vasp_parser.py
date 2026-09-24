"""Real ASE XML observations have identical standalone and workflow conventions."""
from pathlib import Path
import shutil
import ase
import numpy as np
import pytest
from ase.units import GPa
from pymatgen.io.vasp.outputs import Vasprun
from src.utils.dft.vasp_parser import VASPParser


@pytest.mark.parametrize("filename,energy,stress", [
    ("vasprun_pstress.xml", -20.24058451, -165.9956258),
    ("vasprun_dfpt.xml", -6.74587304, 36.08907525),
])
def test_standalone_observations(tmp_path, filename, energy, stress):
    source = Path(ase.__file__).parent / "test/testdata/vasp" / filename
    shutil.copyfile(source, tmp_path / "vasprun.xml")
    parser = VASPParser(tmp_path)
    result = parser.parse_vasprun()
    assert result["final_energy"] == pytest.approx(energy, abs=1e-9)
    np.testing.assert_allclose(result["stress"], np.eye(3) * (-stress * GPa / 10))
    raw = Vasprun(source, parse_potcar_file=False)
    np.testing.assert_array_equal(result["forces"], raw.ionic_steps[-1]["forces"])
    converted = parser.convert_to_matgl_format(result)
    assert converted[0]["energy"] == result["final_energy"]
    for entry, step in zip(converted[1:], raw.ionic_steps[:-1], strict=True):
        np.testing.assert_allclose(entry["structure"].positions, step["structure"].cart_coords)
        assert entry["energy"] == step["e_wo_entrp"]


def test_named_nested_batch_keeps_distinct_source_identities(tmp_path):
    for group, source in [("a", "vasprun_pstress.xml"), ("b", "vasprun_dfpt.xml")]:
        folder = tmp_path / group / "Si.cif"
        folder.mkdir(parents=True)
        shutil.copyfile(Path(ase.__file__).parent / "test/testdata/vasp" / source, folder / "vasprun.xml")
    results = VASPParser(tmp_path).parse_all()
    assert [row["structure_id"] for row in results] == ["a/Si.cif", "b/Si.cif"]
    assert [row["final_energy"] for row in results] == pytest.approx([-20.24058451, -6.74587304])
