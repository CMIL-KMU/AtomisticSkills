"""Explicit v1 multi-geometry, mode and MD contracts for the HTVS registry."""

from copy import deepcopy
from pathlib import Path
import json
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar, Incar


def require(condition, message):
    if not condition:
        raise ValueError(message)


def vector(value, count, label):
    array = np.asarray(value)
    require(array.dtype.kind in "fi", label + " requires numeric values")
    require(
        array.shape == (count, 3) and np.isfinite(array).all(),
        label + " requires finite N by 3 Cartesian values",
    )
    return array


def mobile_projection(values, structure):
    """Orthogonal Cartesian projection onto each atom's mobile lattice span."""
    masks = structure.site_properties.get(
        "selective_dynamics", [[True] * 3] * len(structure)
    )
    projected = []
    for value, mask in zip(values, masks):
        basis = structure.lattice.matrix[np.asarray(mask, dtype=bool)]
        projected.append(np.asarray(value) @ np.linalg.pinv(basis) @ basis)
    return np.asarray(projected)


def attach(prepared, kind, contract):
    """Resolve explicit choices; legacy preparation remains available for comparison."""
    c = deepcopy(contract)
    require(
        c.get("schema") == "vasp-dynamics/v1",
        "Explicit vasp-dynamics/v1 contract required",
    )
    s = Structure.from_dict(prepared["input_structure"])
    order = prepared["atom_order"]
    masks = np.asarray(
        s.site_properties.get("selective_dynamics", [[True] * 3] * len(s))
    )
    require(
        masks.shape == (len(s), 3) and masks.dtype == bool,
        "Boolean atom constraint matrix required",
    )
    common = {"schema", "bindings", "sources"}
    fields = {
        "neb": {
            "intermediate_images",
            "spring",
            "climb",
            "force_tolerance",
            "max_steps",
        },
        "bomd": {
            "ensemble",
            "timestep_fs",
            "equilibration_steps",
            "production_steps",
            "frame_stride",
            "velocities",
            "velocity_unit",
            "temperature_K",
            "smass",
            "langevin_gamma",
            "langevin_gamma_l",
            "pmass",
            "pressure_kbar",
            "random_seed",
        },
        "dimer": {
            "mode",
            "mode_basis",
            "force_tolerance",
            "max_steps",
            "electronic_algorithm",
            "ddr",
            "drotmax",
            "dfnmin",
            "dfnmax",
        },
        "irc": {
            "mode",
            "mode_basis",
            "direction",
            "max_steps",
            "stop",
            "vnorm0",
            "delta0",
            "minstep",
            "maxstep",
        },
    }
    require(
        all(set(b) == {"structure"} for b in c.get("bindings", [])),
        "Invalid ordered geometry input",
    )
    require(
        kind in fields and not set(c) - common - fields[kind], "Unknown dynamics option"
    )
    if kind != "neb":
        require(not c.get("bindings"), "Only NEB accepts additional geometries")
    controlled = {
        "nsteps",
        "nimages",
        "spring",
        "climb",
        "fix_cell",
        "isif",
        "ibrion",
        "ediffg",
        "irc_direction",
        "irc_stop",
        "irc_minstep",
        "timestep",
        "mdalgo",
        "temperature",
        "smass",
        "pmass",
        "ddr",
        "drotmax",
        "dfnmin",
        "dfnmax",
    }
    require(
        not controlled & set(prepared["overrides"]),
        "Dynamics controls belong only in the versioned dynamics contract",
    )
    stage = prepared["stages"][0]
    inc = stage["incar"]
    if kind in ("neb", "dimer", "irc"):
        require(
            type(c.get("max_steps")) is int and c["max_steps"] > 1,
            "Explicit max_steps > 1 required",
        )
        inc.update(NSW=c["max_steps"], ISIF=2)
    if kind in ("neb", "dimer"):
        require(
            type(c.get("force_tolerance")) in (int, float) and c["force_tolerance"] > 0,
            "Positive force_tolerance in eV/angstrom required",
        )
        inc["EDIFFG"] = -c["force_tolerance"]
    if kind == "neb":
        n = c.get("intermediate_images")
        require(
            type(n) is int and 1 <= n <= 97,
            "intermediate_images must be 1..97; total images = IMAGES + 2",
        )
        bindings = c.get("bindings", [])
        require(
            len(bindings) == n + 1,
            "NEB needs all supplied intermediate images and final endpoint",
        )
        images = [s, *[Structure.from_dict(b["structure"]) for b in bindings]]
        for image in images:
            require(
                image.species == s.species
                and np.allclose(
                    image.lattice.matrix, s.lattice.matrix, atol=1e-9, rtol=0
                ),
                "NEB endpoints/images require identical ordered species and cells",
            )
            require(
                image.site_properties == s.site_properties,
                "NEB image constraints/site metadata differ",
            )
        require(
            not np.allclose(images[0].frac_coords, images[-1].frac_coords),
            "NEB endpoints must differ",
        )
        require(
            type(c.get("climb")) is bool
            and type(c.get("spring")) in (int, float)
            and c["spring"] < 0,
            "Explicit climb boolean and negative NEB spring required",
        )
        for a, b in zip(images, images[1:]):
            delta = a.frac_coords - b.frac_coords
            require(
                not np.allclose(delta - np.round(delta), 0),
                "Adjacent NEB images must differ under periodic boundaries",
            )
        inc.update(IMAGES=n, SPRING=c["spring"], LCLIMB=c["climb"])
        stage["images"] = [
            Structure.from_sites([image[i] for i in order]).as_dict()
            for image in images
        ]
        endpoints = []
        for index in (0, len(images) - 1):
            endpoint = deepcopy(stage)
            endpoint.pop("images")
            endpoint.update(
                name="endpoint-" + str(index),
                structure=stage["images"][index],
                input_index=index,
            )
            endpoint["incar"] = {
                k: v for k, v in inc.items() if k not in ("IMAGES", "SPRING", "LCLIMB")
            }
            endpoint["incar"].update(NSW=0, IBRION=-1)
            endpoints.append(endpoint)
        stage["structure"] = stage["images"][1]
        prepared["stages"] = [*endpoints, stage]
    if kind in ("dimer", "irc"):
        require(
            c.get("mode_basis") == "Cartesian-unit-vector",
            "Mode basis must be Cartesian-unit-vector in original atom order",
        )
        mode = vector(c.get("mode"), len(s), "Mode")
        require(
            np.isclose(np.linalg.norm(mode), 1)
            and np.allclose(mode, mobile_projection(mode, s), atol=1e-12, rtol=0),
            "Mode must be normalized and zero on constrained components",
        )
        stage["mode"] = mode[order].tolist()
        if kind == "dimer":
            require(
                c.get("electronic_algorithm") in ("Normal", "Fast"),
                "Dimer requires explicit electronic_algorithm Normal or Fast; legacy ALGO=FALSE cannot establish electronic convergence",
            )
            inc["ALGO"] = c["electronic_algorithm"]
            for key in ("ddr", "dfnmin", "dfnmax"):
                require(
                    type(c.get(key)) in (int, float) and c[key] > 0,
                    "Explicit positive dimer rotation parameters required",
                )
            require(
                type(c.get("drotmax")) is int
                and c["drotmax"] > 0
                and c["dfnmin"] <= c["dfnmax"],
                "Invalid dimer rotation bounds",
            )
            inc.update(
                IBRION=3,
                ICHAIN=2,
                IOPT=2,
                POTIM=0,
                DDR=c["ddr"],
                DROTMAX=c["drotmax"],
                DFNMIN=c["dfnmin"],
                DFNMAX=c["dfnmax"],
            )
        else:
            require(
                type(c.get("direction")) is int and c["direction"] in (-1, 1),
                "IRC direction must be -1 or +1",
            )
            require(
                type(c.get("stop")) is int and c["stop"] > 1,
                "IRC stop must be an integer > 1",
            )
            for key in ("vnorm0", "delta0", "minstep", "maxstep"):
                require(
                    type(c.get(key)) in (int, float) and c[key] > 0,
                    "Explicit positive IRC integration controls required",
                )
            require(c["minstep"] <= c["maxstep"], "IRC time step bounds reversed")
            inc.update(
                IRC_DIRECTION=c["direction"],
                IRC_STOP=c["stop"],
                IRC_VNORM0=c["vnorm0"],
                IRC_DELTA0=c["delta0"],
                IRC_MINSTEP=c["minstep"],
                IRC_MAXSTEP=c["maxstep"],
            )
    if kind == "bomd":
        used = {
            "NVE": set(),
            "NVT": {"temperature_K", "smass"},
            "NPT": {
                "temperature_K",
                "langevin_gamma",
                "langevin_gamma_l",
                "pmass",
                "pressure_kbar",
                "random_seed",
            },
        }
        optional = set.union(*used.values())
        require(
            not (set(c) & optional) - used.get(c.get("ensemble"), set()),
            "Unused ensemble control",
        )
        require(
            c.get("ensemble") in ("NVE", "NVT", "NPT"),
            "Explicit NVE, NVT or NPT ensemble required",
        )
        for key in ("equilibration_steps", "production_steps", "frame_stride"):
            require(
                type(c.get(key)) is int
                and c[key] >= (0 if key == "equilibration_steps" else 1),
                "Explicit MD step counts/stride required",
            )
        require(
            c["frame_stride"] == 1,
            "Dynamics v1 indexes every XML ionic sample; frame_stride must be 1",
        )
        require(
            type(c.get("timestep_fs")) in (int, float) and c["timestep_fs"] > 0,
            "Positive MD timestep_fs required",
        )
        require(
            c.get("velocity_unit") == "angstrom/fs",
            "Initial velocities must use angstrom/fs",
        )
        velocity = vector(c.get("velocities"), len(s), "Velocities")
        require(
            np.allclose(velocity, mobile_projection(velocity, s), atol=1e-12, rtol=0),
            "Constrained lattice-direction velocities must be zero",
        )
        stage["velocities"] = velocity[order].tolist()
        inc.update(
            NSW=c["equilibration_steps"] + c["production_steps"],
            POTIM=c["timestep_fs"],
            NBLOCK=1,
            KBLOCK=1,
            ISIF=2,
            MDALGO=2,
            SMASS=-3,
        )
        for key in ("TEBEG", "TEEND", "PMASS", "LANGEVIN_GAMMA", "LANGEVIN_GAMMA_L"):
            inc.pop(key, None)
        if c["ensemble"] != "NVE":
            require(
                type(c.get("temperature_K")) in (int, float) and c["temperature_K"] > 0,
                "Positive temperature_K required",
            )
            inc.update(TEBEG=c["temperature_K"], TEEND=c["temperature_K"])
        if c["ensemble"] == "NVT":
            require(
                type(c.get("smass")) in (int, float) and c["smass"] >= 0,
                "NVT requires explicit nonnegative Nose smass",
            )
            inc["SMASS"] = c["smass"]
        if c["ensemble"] == "NPT":
            for key in ("pmass", "langevin_gamma_l"):
                require(
                    type(c.get(key)) in (int, float) and c[key] > 0,
                    "NPT requires explicit positive lattice damping and mass",
                )
            gamma = c.get("langevin_gamma", [])
            require(
                len(gamma) == len(prepared["potcar_symbols"])
                and all(type(x) in (int, float) and x > 0 for x in gamma),
                "Langevin damping must follow PAW species order",
            )
            require(
                type(c.get("pressure_kbar")) in (int, float),
                "NPT requires pressure_kbar",
            )
            seed = c.get("random_seed", [])
            require(
                len(seed) == 3 and all(type(x) is int and x >= 0 for x in seed),
                "NPT requires explicit three-integer RANDOM_SEED",
            )
            inc.pop("SMASS")
            inc.update(
                ISIF=3,
                MDALGO=3,
                LANGEVIN_GAMMA=gamma,
                LANGEVIN_GAMMA_L=c["langevin_gamma_l"],
                PMASS=c["pmass"],
                PSTRESS=c["pressure_kbar"],
                RANDOM_SEED=Incar.from_str("RANDOM_SEED = "+" ".join(map(str,seed)))["RANDOM_SEED"],
            )
    stage["dynamics"] = dict(c, kind=kind)
    prepared["dynamics"] = c
    return prepared


def setup(folder, stage):
    """Additional native inputs; never randomize the supplied coordinates."""
    if "images" in stage:
        for i, image in enumerate(stage["images"]):
            directory = folder / f"{i:02d}"
            directory.mkdir()
            Poscar(Structure.from_dict(image)).write_file(directory / "POSCAR")
            for name in ("INCAR", "KPOINTS"):
                (directory / name).write_bytes((folder / name).read_bytes())
    if "velocities" in stage:
        poscar = Poscar.from_file(folder / "POSCAR", read_velocities=False)
        Poscar(poscar.structure, velocities=stage["velocities"]).write_file(
            folder / "POSCAR"
        )
    if "mode" in stage:
        text = "\n".join(" ".join(map(str, row)) for row in stage["mode"]) + "\n"
        if stage["dynamics"]["kind"] == "dimer":
            (folder / "MODECAR").write_text(text)
        else:
            with (folder / "POSCAR").open("a") as stream:
                stream.write(
                    "! unstable direction optimized by the dimer method\n" + text
                )


def neb_observables(path):
    from pymatgen.io.vasp.outputs import Outcar

    outcar = Outcar(path)
    outcar.read_neb()
    return dict(
        sigma_zero_energy_eV=outcar.data["energy"],
        tangent_force_eV_per_angstrom=outcar.data.get("tangent_force"),
    )


def parse(directory, stage):
    """Delegate numerical parsing, adding family-specific completion evidence."""
    from src.simulations.vasp.parser import parse_stage

    directory = Path(directory)
    if "images" in stage:
        reports = []
        for i in range(1, len(stage["images"]) - 1):
            member = directory / f"{i:02d}"
            report = parse_stage(member)
            report["neb"] = neb_observables(member / "OUTCAR")
            require(
                np.isclose(
                    report["neb"]["sigma_zero_energy_eV"],
                    report["steps"][-1]["energies"]["e_0_energy"],
                    atol=1e-5,
                ),
                "NEB XML and OUTCAR final energies disagree",
            )
            text = (member / "OUTCAR").read_text()
            report["converged"] = (
                report["electronic_converged"]
                and "reached required accuracy" in text
                and len(report["steps"]) < stage["incar"]["NSW"]
            )
            reports.append(report)
        return dict(
            images=reports,
            converged=all(r["converged"] for r in reports),
            steps=reports[-1]["steps"],
        )
    report = parse_stage(directory)
    c = stage.get("dynamics")
    if not c:
        return report
    steps = report["steps"]
    kind = c["kind"]
    if kind == "bomd":
        report["converged"] = (
            report["electronic_converged"] and len(steps) == stage["incar"]["NSW"]
        )
    elif kind == "dimer":
        path = directory / "DIMCAR"
        require(
            path.read_text().splitlines()[0].split()
            == ["Step", "Force", "Torque", "Energy", "Curvature", "Angle"],
            "Unsupported VTST DIMCAR header; verify selected patch output layout",
        )
        data = np.genfromtxt(
            path, skip_header=1, missing_values="---", filling_values=np.nan, ndmin=2
        )
        require(
            data.shape[1] == 6 and len(data) and np.isfinite(data[:, [0, 1, 3]]).all(),
            "VTST DIMCAR requires finite step, force and energy",
        )
        rotation = len(data) - 1
        if not np.isfinite(data[-1]).all():
            require(
                len(data) > 1
                and np.isnan(data[-1, [2, 4, 5]]).all()
                and data[-1, 0] - data[-2, 0] in (0, 1)
                and np.isclose(data[-1, 3], data[-2, 3], atol=1e-5),
                "Invalid VTST terminal placeholder row",
            )
            rotation -= 1
        require(
            np.isfinite(data[: rotation + 1]).all(),
            "Missing rotational DIMCAR evidence",
        )
        residual = mobile_projection(
            steps[-1]["forces"], Structure.from_dict(report["input_structure"])
        )
        report["dimer"] = dict(
            curvature=float(data[rotation, 4]),
            curvature_step=int(data[rotation, 0]),
            force=float(data[-1, 1]),
            classification="converged saddle candidate; Hessian index unverified",
        )
        report["converged"] = (
            report["electronic_converged"]
            and report["ionic_converged"]
            and data[rotation, 4] < 0
            and "reached required accuracy" in (directory / "OUTCAR").read_text()
            and np.max(np.linalg.norm(residual, axis=1)) <= c["force_tolerance"]
        )
    elif kind == "irc":
        values = [s["energies"]["e_fr_energy"] for s in steps]
        stop = c["stop"]
        complete = (
            len(values) > stop
            and all(a < b for a, b in zip(values[-stop - 1 : -1], values[-stop:]))
            and len(values) < c["max_steps"]
        )
        report["converged"] = report["electronic_converged"] and complete
        report["irc"] = dict(
            direction=c["direction"],
            endpoint_index=int(np.argmin(values)),
            endpoint="lowest free-energy sampled point; uphill tail retained",
            scope="one saddle-to-basin branch, not a two-minimum reaction barrier",
        )
    report["converged"] = bool(report["converged"])
    return json.loads(json.dumps(report, allow_nan=False))
