"""Versioned HTVS scientific declarations; no HTVS runtime or site configuration."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
import os
import yaml
from importlib.metadata import version
from importlib.util import find_spec
from jinja2 import StrictUndefined
from jinja2.sandbox import SandboxedEnvironment

ROOT = Path(os.environ["ATOMISTIC_VASP_PROTOCOLS"]) if "ATOMISTIC_VASP_PROTOCOLS" in os.environ else None


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def kpoint_data(value):
    data = value.as_dict()
    if data["labels"]:
        data["labels"] = [label or "" for label in data["labels"]]
    return data


def implementation():
    """Measure source and scientific dependencies independently of any consumer."""
    from src.simulations.identity import identity
    return identity("vasp")


def registry():
    if ROOT is None:
        raise ValueError("An explicit protocol directory is required")
    data = yaml.safe_load((ROOT / "registry.yaml").read_text())
    aliases = {}
    for profile in data["profiles"]:
        for name in [profile["id"], *profile["aliases"]]:
            if name in aliases:
                raise ValueError("Duplicate protocol identity or alias: " + name)
            aliases[name] = profile
    return data, aliases


def resolve(name):
    data, aliases = registry()
    if name not in aliases:
        raise ValueError("Unknown VASP protocol")
    return data, deepcopy(aliases[name])


def prepare(name, structure, overrides, *, selections, dynamics=None, ldau_preset=None):
    """Freeze scientific stages for independent calculation (no POTCAR required)."""
    from pymatgen.core import Structure, Species
    from pymatgen.io.vasp.inputs import Incar, Kpoints

    data, profile = resolve(name)
    if set(selections) != set(profile["ambiguities"]) or any(
        v != "htvs-effective" for v in selections.values()
    ):
        raise ValueError(
            "Explicit htvs-effective selection required for: "
            + ", ".join(profile["ambiguities"])
        )
    text = "".join(
        data["fragments"][k]
        for stage in profile["stages"]
        for k in data["templates"][stage["template"]]
    )
    optional = set(re.findall(r"details.get\(['\"]([a-z_]+)", text))
    allowed = (
        set(re.findall(r"details\.([a-z_]+)", text))
        | optional
        | {
            "magmoms",
            "oxidation_states",
            "surface_atoms",
            "adsorbate_atoms",
            "irc_direction",
            "kpoints",
            "kppa",
            "kppvol",
        }
    )
    allowed -= {"get", "magmom", "ldauu", "ldauj", "ldaul", "lmaxmix", "nelect"}
    if ldau_preset is not None:
        allowed -= {"ldau", "ldautype", "ldauprint"}
    allowed.add("fix_cell")
    if profile["kind"] == "dielectric":
        allowed.add("ismear")
    if "_surf_" not in profile["id"] or profile["kind"] != "opt":
        allowed -= {"surface_atoms", "adsorbate_atoms"}
    if profile["kind"] != "irc":
        allowed.discard("irc_direction")
    if set(overrides) - allowed or set(overrides) & set(data["machine_keys"]):
        raise ValueError("Unknown or machine-owned scientific override")
    details = dict(profile["defaults"], **overrides)
    if dynamics is not None and profile["kind"] == "irc":
        details["irc_direction"] = dynamics.get("direction")
    if profile["kind"] == "dielectric" and not {"fix_cell", "ismear"} <= set(overrides):
        raise ValueError(
            "Dielectric helper defaults are absent in HTVS; explicitly select fix_cell and ismear"
        )
    selected_mesh = set(overrides) & {"kpoints", "kppa", "kppvol"}
    if len(selected_mesh) > 1:
        raise ValueError("Select one k-point convention")
    if selected_mesh:
        for key in {"kpoints", "kppa", "kppvol"} - selected_mesh:
            details.pop(key, None)
    for key, value in overrides.items():
        default = profile["defaults"].get(key)
        if (
            default is not None
            and type(value) is not type(default)
            and not (type(default) is float and type(value) is int)
        ):
            raise ValueError("Override type mismatch: " + key)
    if details.get("bader") or details.get("charge") or details.get("nelect"):
        raise ValueError(
            "Charged/Bader protocols require a separately qualified contract"
        )
    original = Structure.from_dict(structure)
    if not all(original.lattice.pbc):
        raise ValueError("VASP requires full periodicity; use an explicit periodic vacuum cell")
    s = original.copy()
    if not s.is_ordered or not len(s):
        raise ValueError("Ordered periodic structure required")
    order = sorted(range(len(s)), key=lambda i: s[i].specie.symbol)
    s = Structure.from_sites([s[i] for i in order])
    elements = list(dict.fromkeys(site.specie.symbol for site in s))
    if any(e not in data["tables"]["potcar"] for e in elements):
        raise ValueError("No HTVS PAW selection for element")
    if details.get("kpoints"):
        mesh = details["kpoints"]
        if len(mesh) != 3 or any(type(x) is not int or x < 1 for x in mesh):
            raise ValueError("Three positive mesh integers required")
        kpoints = Kpoints.gamma_automatic(mesh)
    elif details.get("kppa"):
        kpoints = Kpoints.automatic_density(s, details["kppa"])
    else:
        kpoints = Kpoints.automatic_density_by_vol(s, details.get("kppvol", 64))
    details["ismear"] = 0 if min(kpoints.kpts[0]) < 4 else details.get("ismear", 0)
    details["isif"] = details.get("isif", 2 if details.get("fix_cell", False) else 3)
    details["nelect"] = None
    moments = details.get("magmoms", original.site_properties.get("magmom"))
    oxidation = details.get(
        "oxidation_states", original.site_properties.get("oxidation_states")
    )
    if oxidation is not None and (
        len(oxidation) != len(s)
        or any(type(v) not in (int, float) for v in oxidation)
        or abs(sum(oxidation)) > 1e-8
    ):
        raise ValueError(
            "Explicit oxidation states must match atom order and neutral charge"
        )
    if moments is not None and len(moments) != len(s):
        raise ValueError("MAGMOM must match input atom ordering")
    details["magmom"] = " ".join(
        str(
            moments[i]
            if moments is not None
            else data["tables"]["magmom"].get(
                str(
                    original[i].specie
                    if oxidation is None
                    else Species(original[i].specie.symbol, oxidation[i])
                ),
                data["tables"]["magmom"].get(original[i].specie.symbol, 0),
            )
        )
        for i in order
    )
    for key, default in [("ldauu", 0), ("ldauj", 0), ("ldaul", -1)]:
        details[key] = " ".join(
            str(data["tables"][key].get(e, default)) for e in elements
        )
    details["lmaxmix"] = (
        6
        if max(site.specie.Z for site in s) > 56
        else 4
        if max(site.specie.Z for site in s) > 20
        else 2
    )
    if "_surf_" in name and profile["kind"] == "opt":
        masks = [details.get(k) for k in ("surface_atoms", "adsorbate_atoms")]
        if any(
            not isinstance(m, list)
            or len(m) != len(s)
            or any(type(v) is not bool for v in m)
            for m in masks
        ):
            raise ValueError(
                "Explicit surface/adsorbate masks in original atom order required"
            )
        s.add_site_property(
            "selective_dynamics", [[masks[0][i] or masks[1][i]] * 3 for i in order]
        )
    env = SandboxedEnvironment(undefined=StrictUndefined)
    stages = []
    for stage in profile["stages"]:
        template = "".join(
            data["fragments"][k] for k in data["templates"][stage["template"]]
        )
        incar = dict(
            Incar.from_str(
                env.from_string(template).render(jobspec={"details": details})
            )
        )
        if set(incar) & set(data["parallel_keys"]):
            raise ValueError("Machine parallelism leaked into scientific protocol")
        stage_kpoints = kpoints
        if (
            stage["name"] == "dos"
            and not details.get("kpoints")
            and not details.get("kppa")
        ):
            stage_kpoints = Kpoints.automatic_density_by_vol(
                s, details.get("kppvol", 64) * 4
            )
        if stage["name"] == "bs":
            from pymatgen.symmetry.bandstructure import HighSymmKpath

            points, labels = HighSymmKpath(s).get_kpoints(
                line_density=details.get("kpath_density", 20),
                coords_are_cartesian=False,
            )
            stage_kpoints = Kpoints(
                comment="HTVS reciprocal path",
                style=Kpoints.supported_modes.Reciprocal,
                num_kpts=len(points),
                kpts=[[round(float(x), 5) for x in point] for point in points],
                kpts_weights=[1.0] * len(points),
                labels=labels,
            )
        stages.append(
            dict(
                name=stage["name"],
                incar=incar,
                kpoints=kpoint_data(Kpoints.from_str(str(stage_kpoints))),
            )
        )
    prepared = dict(
        schema="atomistic-skills.vasp-input/v1",
        protocol=profile["id"],
        registry_sha256=digest(data),
        implementation=implementation(),
        htvs_revision=data["source"]["revision"],
        source=profile["source_hashes"],
        overrides=overrides,
        selections=selections,
        details=details,
        input_structure=structure,
        structure=s.as_dict(),
        atom_order=order,
        potcar_symbols=[data["tables"]["potcar"][e] for e in elements],
        stages=stages,
    )
    if dynamics is not None:
        from src.simulations.vasp.dynamics import attach

        attach(prepared, profile["kind"], dynamics)
    if ldau_preset is not None:
        from src.simulations.vasp.ldau import apply

        apply(prepared, ldau_preset)
    return json.loads(json.dumps(prepared, allow_nan=False))

