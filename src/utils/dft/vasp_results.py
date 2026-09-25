"""VASP observations with explicit conventions and absent optional quantities.

Full stress is required only for ISIF >= 2; trace-only ISIF=1 is not a tensor.
Spectral objects are optional unless a consuming calculation requests them.
Runtime extension parameters retain their OUTCAR origin separately from XML.
"""

from pathlib import Path
from copy import copy
import json
import re
from monty.json import MontyEncoder
import numpy as np
from ase.units import GPa
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints
def kpoint_data(value):
    if value is None:
        return None
    data = value.as_dict()
    if data["labels"]:
        data["labels"] = [label or "" for label in data["labels"]]
    return data


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
    # ISIF=1 provides only a trace, not a usable full stress tensor.
    parameters = dict(v.parameters, **v.incar)
    full_stress = parameters.get("ISIF", 0 if parameters.get("IBRION") == 0 or parameters.get("LHFCALC") else 2) >= 2
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
                **({"stress": (-np.asarray(step["stress"]) * GPa / 10).tolist()} if full_stress else {}),
            )
        )
    if not steps:
        raise ValueError("VASP output contains no ionic observations")
    return steps


def parse_stage(directory):
    input_incar = Incar.from_file(Path(directory) / "INCAR")
    if not (Path(directory) / "KPOINTS").is_file() and (
        type(input_incar.get("KSPACING")) not in (int, float) or input_incar["KSPACING"] <= 0
    ):
        raise ValueError("Output lacks explicit KPOINTS or positive KSPACING")
    v = Vasprun(Path(directory) / "vasprun.xml", parse_potcar_file=False, parse_projected_eigen=False)
    outcar = Outcar(Path(directory) / "OUTCAR")
    steps = ionic_observables(v)
    # Reuse pymatgen's LEPSILON/Exact/ordinary SCF rules for every ionic sample.
    # Response iterations are not additional unconverged SCF iterations.
    sample = copy(v)
    electronic_converged = True
    for step in v.ionic_steps:
        sample.ionic_steps = [step]
        electronic_converged &= bool(sample.converged_electronic)
    result = dict(
        units=UNITS,
        version=v.vasp_version,
        converged=bool(v.converged),
        electronic_converged=electronic_converged,
        ionic_converged=bool(v.converged_ionic),
        initial_structure=v.initial_structure.as_dict(),
        steps=steps,
        incar=dict(v.incar),
        engine_parameters=dict(v.parameters),
        input_incar=dict(input_incar),
        input_kpoints=kpoint_data(Kpoints.from_file(Path(directory) / "KPOINTS"))
        if (Path(directory) / "KPOINTS").is_file() else None,
        input_structure=Poscar.from_file(
            Path(directory) / "POSCAR"
        ).structure.as_dict(),
        potcar_titles=v.potcar_symbols,
        run_stats=outcar.run_stats,
    )
    result["converged"] = result["converged"] and result["electronic_converged"]
    text = (Path(directory) / "OUTCAR").read_text()
    controls = {}
    vtst = re.search(r"^\s*VTST: version\s+([\d.]+)", text, re.MULTILINE)
    if vtst:
        result["extensions"] = {"vtst": vtst[1]}
    if v.parameters.get("IBRION", v.incar.get("IBRION")) == 40:
        # XML omits these in VASP 6.4.1. Read the engine's runtime block,
        # never the echoed user INCAR at the beginning of OUTCAR.
        marker = "DAMPED VELOCITY VERLET ALGORITHM:"
        if marker in text:
            block = text.split(marker, 1)[1].split("IRC (A):", 1)[0]
            controls = {k: value for k, value in Incar.from_str(block).items()
                        if k in {"IRC_DIRECTION", "IRC_STOP", "IRC_MINSTEP", "IRC_MAXSTEP", "IRC_VNORM0", "IRC_DELTA0"}}
    # VTST prints its own runtime controls rather than XML INCAR entries.
    for key, pattern in {
        "ICHAIN": r"^\s*CHAIN: Read ICHAIN\s+(\d+)",
        "DROTMAX": r"^\s*Dimer: RotMax\s+(\d+)",
        "DDR": r"^\s*Dimer:\s+dR\s+([\d.Ee+\-]+)",
        "DFNMAX": r"^\s*Dimer:\s+FNMax\s+([\d.Ee+\-]+)",
        "DFNMIN": r"^\s*Dimer:\s+FNMin\s+([\d.Ee+\-]+)",
    }.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            controls[key] = float(match[1])
    if re.search(r"^\s*OPT: Using Conjugate-Gradient optimizer\s*$", text, re.MULTILINE):
        controls["IOPT"] = 2
    if controls:
        result["outcar_parameters"] = controls
    if getattr(v, "tdos", None) is not None:
        result["dos"] = dict(
            energies=v.tdos.energies.tolist(),
            densities={str(int(k)): x.tolist() for k, x in v.tdos.densities.items()},
            efermi=v.efermi,
        )
    if getattr(v, "eigenvalues", None) is not None:
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
