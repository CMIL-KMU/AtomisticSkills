"""Named scientific Hubbard presets, reusing the single generated HTVS table source."""

import math
import hashlib
from pathlib import Path
import yaml
from src.simulations.vasp import protocols
from src.simulations.vasp.protocols import digest, registry


def apply(prepared, selection):
    from pymatgen.core import Structure

    if not isinstance(selection, dict) or set(selection) != {"name", "version"}:
        raise ValueError("Explicit LDAU preset name and version required")
    catalog = yaml.safe_load((protocols.ROOT / "ldau.yaml").read_text())
    if (
        set(catalog) != {"schema", "presets"}
        or catalog["schema"] != "anvil.ldau-presets/v1"
        or selection["name"] not in catalog["presets"]
    ):
        raise ValueError("Unknown LDAU preset")
    preset = catalog["presets"][selection["name"]]
    required = {"version", "mode", "functionals", "paw_family", "reference_policy"}
    if (
        set(preset)
        - required
        - {"table", "missing", "no_u_elements", "ldautype", "values"}
        or not required <= set(preset)
        or preset["version"] != selection["version"]
    ):
        raise ValueError("Invalid preset fields/version")
    stages = prepared["stages"]
    functionals = {s["incar"].get("METAGGA", "PBE").upper() for s in stages}
    if (
        len(functionals) != 1
        or not functionals <= set(preset["functionals"])
        or any(
            s["incar"].get("LHFCALC") or s["incar"].get("GGA", "PE") != "PE"
            for s in stages
        )
    ):
        raise ValueError("LDAU preset is incompatible with effective functional")
    data = registry()[0]
    tables = data["tables"]
    structure = Structure.from_dict(prepared["structure"])
    elements = list(dict.fromkeys(site.specie.symbol for site in structure))
    mode = preset["mode"]
    if mode not in {"legacy", "disabled", "ueff", "separate"}:
        raise ValueError("Unknown Hubbard U/J semantics")
    if (
        mode == "ueff"
        and preset.get("ldautype") != 2
        or mode == "separate"
        and preset.get("ldautype") not in (1, 2, 4)
    ):
        raise ValueError("Incompatible LDAUTYPE and U/J semantics")
    if (
        preset.get("missing", "reject") not in {"reject", "explicit-legacy-zero"}
        or preset.get("missing") == "explicit-legacy-zero"
        and mode != "legacy"
    ):
        raise ValueError(
            "Missing-element policy must be explicit; legacy fallback is compatibility-only"
        )
    values = preset.get("values", {})
    if preset.get("table") == "htvs":
        if values:
            raise ValueError("Conflicting table and explicit Hubbard values")
        values = {
            e: dict(l=tables["ldaul"][e], u=tables["ldauu"][e], j=tables["ldauj"][e])
            for e in tables["ldauu"]
        }
    elif preset.get("table") is not None:
        raise ValueError("Unknown Hubbard table reference")
    if (
        set(values) & set(preset.get("no_u_elements", []))
        or mode == "disabled"
        and values
    ):
        raise ValueError("Conflicting corrected and explicitly uncorrected species")
    arrays = dict(LDAUL=[], LDAUU=[], LDAUJ=[])
    oxidation = prepared["details"].get(
        "oxidation_states",
        Structure.from_dict(prepared["input_structure"]).site_properties.get(
            "oxidation_states"
        ),
    )
    original = Structure.from_dict(prepared["input_structure"])
    for element in elements:
        entry = values.get(element)
        if isinstance(entry, dict) and "oxidation" in entry:
            states = {
                (
                    oxidation[i]
                    if oxidation is not None
                    else getattr(site.specie, "oxi_state", None)
                )
                for i, site in enumerate(original)
                if site.specie.symbol == element
            }
            states = {format(v, "g") if v is not None else "undefined" for v in states}
            if (
                set(entry) != {"oxidation"}
                or len(states) != 1
                or next(iter(states)) not in entry["oxidation"]
            ):
                raise ValueError(
                    "Explicit unambiguous oxidation mapping required per POTCAR species"
                )
            entry = entry["oxidation"][next(iter(states))]
        if (
            mode == "disabled"
            or entry is None
            and (
                preset.get("missing") == "explicit-legacy-zero"
                or element in preset.get("no_u_elements", [])
            )
        ):
            entry = dict(l=-1, u=0, j=0)
        if entry is None or set(entry) != {"l", "u", "j"}:
            raise ValueError(
                "Undefined species or ambiguous U_eff versus U/J mapping: " + element
            )
        if (
            not isinstance(entry["l"], int)
            or isinstance(entry["l"], bool)
            or entry["l"] not in (-1, 0, 1, 2, 3)
            or any(
                type(entry[k]) not in (int, float)
                or not math.isfinite(entry[k])
                or entry[k] < 0
                for k in ("u", "j")
            )
        ):
            raise ValueError("Invalid Hubbard l/U/J values")
        if (
            mode == "ueff"
            and entry["j"] != 0
            or entry["l"] == -1
            and (entry["u"] or entry["j"])
        ):
            raise ValueError("Contradictory U_eff/J or inactive orbital correction")
        for key, field in [("LDAUL", "l"), ("LDAUU", "u"), ("LDAUJ", "j")]:
            arrays[key].append(entry[field])
    for stage in stages:
        incar = stage["incar"]
        if mode == "disabled":
            for key in list(incar):
                if key.startswith("LDAU"):
                    del incar[key]
            incar["LDAU"] = False
        elif mode != "legacy":
            incar.update(arrays, LDAU=True, LDAUTYPE=preset["ldautype"])
    for key, value in arrays.items():
        prepared["details"][key.lower()] = " ".join(str(v) for v in value)
    prepared["schema"] = "atomistic-skills.vasp-input/v2"
    prepared["ldau_preset"] = dict(
        selection=selection,
        definition=preset,
        sha256=digest(dict(preset=preset, tables=values)),
        catalog_sha256=digest(catalog),
        source={k: v for k, v in data["source"].items() if k != "files"},
        elements=elements,
        arrays=arrays,
        effective=[
            {k: v for k, v in s["incar"].items() if k.startswith("LDAU")}
            for s in stages
        ],
    )
    prepared["implementation"]["files"]["ldau.py"] = hashlib.sha256(
        Path(__file__).read_bytes()
    ).hexdigest()
    return prepared
