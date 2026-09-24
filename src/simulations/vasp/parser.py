"""VASP numerical extraction with explicit energy and stress conventions."""

from pathlib import Path
import json
from monty.json import MontyEncoder
import numpy as np
from ase.units import GPa
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints
from src.simulations.vasp.protocols import kpoint_data

UNITS = dict(
    energy="eV",
    forces="eV/angstrom",
    stress="eV/angstrom3",
    stress_sign="tensile",
    energy_observable="ionic_steps.e_wo_entrp",
)


def ionic_observables(v):
    """One energy/stress convention for standalone parsing and workflow consumers."""
    steps = []
    for step in v.ionic_steps:
        steps.append(
            dict(
                structure=step["structure"].as_dict(),
                energy=step["e_wo_entrp"],
                energies={
                    k: step[k]
                    for k in ("e_wo_entrp", "e_fr_energy", "e_0_energy")
                    if k in step
                },
                forces=np.asarray(step["forces"]).tolist(),
                stress=(-np.asarray(step["stress"]) * GPa / 10).tolist(),
            )
        )
    if not steps:
        raise ValueError("VASP output contains no ionic observations")
    return steps


def parse_stage(directory):
    v = Vasprun(Path(directory) / "vasprun.xml", parse_potcar_file=False, parse_projected_eigen=False)
    outcar = Outcar(Path(directory) / "OUTCAR")
    steps = ionic_observables(v)
    result = dict(
        units=UNITS,
        version=v.vasp_version,
        converged=bool(v.converged),
        electronic_converged=bool(v.converged_electronic)
        and all(
            len(step.get("electronic_steps", [])) < v.parameters.get("NELM", 60)
            for step in v.ionic_steps
        ),
        ionic_converged=bool(v.converged_ionic),
        initial_structure=v.initial_structure.as_dict(),
        steps=steps,
        incar=dict(v.incar),
        engine_parameters=dict(v.parameters),
        input_incar=dict(Incar.from_file(Path(directory) / "INCAR")),
        input_kpoints=kpoint_data(Kpoints.from_file(Path(directory) / "KPOINTS")),
        input_structure=Poscar.from_file(
            Path(directory) / "POSCAR"
        ).structure.as_dict(),
        potcar_titles=v.potcar_symbols,
        run_stats=outcar.run_stats,
    )
    if v.tdos is not None:
        result["dos"] = dict(
            energies=v.tdos.energies.tolist(),
            densities={str(int(k)): x.tolist() for k, x in v.tdos.densities.items()},
            efermi=v.efermi,
        )
    if v.eigenvalues is not None:
        result["bands"] = dict(
            kpoints=v.actual_kpoints,
            values={str(int(k)): x.tolist() for k, x in v.eigenvalues.items()},
            efermi=v.efermi,
        )
    if v.incar.get("LEPSILON"):
        result["dielectric"] = {
            key: np.asarray(getattr(v, key)).tolist()
            for key in ("epsilon_static", "epsilon_static_wolfe", "epsilon_ionic")
            if len(getattr(v, key))
        }
        if "epsilon_static" not in result["dielectric"]:
            raise ValueError("Missing requested dielectric tensor")
    if v.incar.get("IBRION") in (5, 6) or hasattr(v, "force_constants"):
        from ase.calculators.vasp import Vasp
        from pymatgen.io.ase import AseAtomsAdaptor
        import xml.etree.ElementTree as ET

        atoms = AseAtomsAdaptor.get_atoms(
            Poscar.from_file(Path(directory) / "POSCAR").structure
        )
        types = ET.parse(
            Path(directory) / "vasprun.xml"
        ).findall(
            "./atominfo/array[@name='atomtypes']/set/rc"  # XML selector, not a filesystem/site path
        )
        masses = [float(row[2].text) for row in types for _ in range(int(row[0].text))]
        if len(masses) != len(atoms) or atoms.constraints:
            raise ValueError(
                "Frequency v1 requires explicit PAW masses and all atoms free"
            )
        atoms.set_masses(masses)
        calc = Vasp(directory=str(directory))
        calc.atoms = atoms
        calc.sort = list(range(len(atoms)))
        vibrations = calc.get_vibrations()
        frequencies = vibrations.get_frequencies()
        result["vibrations"] = dict(
            hessian=vibrations.get_hessian_2d().tolist(),
            frequencies=[
                float(f.real) if abs(f.imag) < 1e-8 else -abs(float(f.imag))
                for f in frequencies
            ],
            masses=masses,
            hessian_unit="eV/angstrom2",
            frequency_unit="1/cm",
            imaginary="negative frequency",
        )
    return json.loads(json.dumps(result, cls=MontyEncoder, allow_nan=False))
