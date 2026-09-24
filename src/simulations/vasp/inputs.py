"""Materialize scientific inputs from explicit, caller-supplied PAW datasets."""
from pathlib import Path
import shutil
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from src.simulations.vasp.dynamics import setup
from src.simulations.vasp.protocols import kpoint_data


def write_inputs(folder, structure, stage, potcar_paths, parallel=None, previous=None):
    """Write inputs into an empty directory; scheduling and binary choice are external."""
    folder = Path(folder)
    if folder.exists() and any(folder.iterdir()):
        raise FileExistsError("VASP input directory must be empty")
    folder.mkdir(parents=True, exist_ok=True)
    potcar = Potcar()
    for path in potcar_paths:
        potcar.extend(Potcar.from_file(path))
    incar = Incar(dict(stage["incar"], **(parallel or {})))
    kpoints = Kpoints.from_dict(stage["kpoints"])
    VaspInput(incar, kpoints, Poscar(Structure.from_dict(structure)), potcar).write_input(folder)
    setup(folder, stage)
    if previous is not None and incar.get("ICHARG") in (1, 11):
        shutil.copyfile(Path(previous) / "CHGCAR", folder / "CHGCAR")
    return dict(incar=dict(Incar.from_file(folder / "INCAR")), kpoints=kpoint_data(kpoints),
                structure=Poscar.from_file(folder / "POSCAR").structure.as_dict())
