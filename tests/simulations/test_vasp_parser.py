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


@pytest.mark.parametrize("isif", [0, 1, 2])
def test_stress_requires_full_tensor_calculation(isif):
    from types import SimpleNamespace
    from pymatgen.core import Structure, Lattice
    from src.utils.dft.vasp_results import ionic_observables
    step = dict(structure=Structure(Lattice.cubic(3), ["Si"], [[0, 0, 0]]),
                e_wo_entrp=-1., forces=[[0., 0., 0.]])
    v = SimpleNamespace(incar={"ISIF": isif}, parameters={}, ionic_steps=[step])
    if isif == 2:
        with pytest.raises(KeyError, match="stress"):
            ionic_observables(v)
    else:
        assert "stress" not in ionic_observables(v)[0]
    step["stress"] = np.eye(3)
    result = ionic_observables(v)[0]
    assert ("stress" in result) is (isif == 2)


def test_missing_stress_standalone_conversion(tmp_path):
    import xml.etree.ElementTree as ET
    source = Path(ase.__file__).parent / "test/testdata/vasp/vasprun_pstress.xml"
    tree = ET.parse(source)
    for elem in tree.findall(".//i[@name='ISIF']"):
        elem.text = "0"
    for calc in tree.findall(".//calculation"):
        stress = calc.find("varray[@name='stress']")
        if stress is not None:
            calc.remove(stress)
    tree.write(tmp_path / "vasprun.xml")
    parser = VASPParser(tmp_path)
    result = parser.parse_vasprun()
    assert result["stress"] is None
    assert all(row["stress"] is None for row in parser.convert_to_matgl_format(result))


def test_optional_spectral_objects_and_runtime_controls(tmp_path, monkeypatch):
    from pymatgen.io.vasp.inputs import Poscar
    from src.utils.dft import vasp_results
    root = Path(ase.__file__).parent / 'test/testdata/vasp'
    v = Vasprun(root / 'vasprun_pstress.xml', parse_potcar_file=False)
    v.incar.write_file(tmp_path / 'INCAR')
    v.kpoints.write_file(tmp_path / 'KPOINTS')
    Poscar(v.initial_structure).write_file(tmp_path / 'POSCAR')
    for name in ('tdos', 'eigenvalues'):
        if hasattr(v, name):
            delattr(v, name)
    v.parameters['IBRION'] = 40
    shutil.copyfile(root / 'OUTCAR_example_1', tmp_path / 'OUTCAR')
    with (tmp_path / 'OUTCAR').open('a') as out:
        out.write('\n IRC_STOP = 99\n DAMPED VELOCITY VERLET ALGORITHM:\n Specific input parameters:\n'
                  ' IRC_STOP = 3\n IRC_DIRECTION = 1\n IRC_MINSTEP = 0.0500\n IRC_MAXSTEP = 3.0000\n'
                  ' IRC_VNORM0 = 0.0030\n IRC_DELTA0 = 0.0015\n IRC (A): 0.000 E(eV): -1.0\n')
    monkeypatch.setattr(vasp_results, 'Vasprun', lambda *args, **kwargs: v)
    report = vasp_results.parse_stage(tmp_path)
    assert 'dos' not in report and 'bands' not in report
    assert report['outcar_parameters']['IRC_STOP'] == 3
    assert report['outcar_parameters']['IRC_DELTA0'] == .0015
    assert 'IRC_STOP' not in report['incar']


def test_fixed_charge_cannot_become_training_labels(tmp_path):
    with pytest.raises(ValueError, match="Fixed-charge"):
        VASPParser(tmp_path).convert_to_matgl_format({"incar": {"ICHARG": 11}})
