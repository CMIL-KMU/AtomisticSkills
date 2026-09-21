"""Single Anvil material policy and versioned contracts; upstream tools are intact."""

import hashlib
import json
import os
from pathlib import Path
from atomistic_analysis import tools

CATALOG = json.loads(Path(__file__).with_name("capabilities.json").read_text())
OPERATIONS = tuple(CATALOG["approved_operations"])
SCHEMA = "anvil.material-request/v2"

# Root input ports accepting registered scientific data. Scalar ports also accept
# individual list indices; arbitrary object traversal and code dispatch are absent.
PORTS = {
    "xrd": {"structure": "structure"},
    "supercell": {"structure": "structure"},
    "substitute": {"structure": "structure"},
    "slab": {"structure": "structure"},
    "symmetry": {"structure": "structure"},
    "structure-match": {"first": "structure", "second": "structure"},
    "rdf": {"frames": "frames"},
    "probability-density": {"frames": "frames"},
    "ideal-gas": {
        "symbols": "symbols",
        "positions_angstrom": "positions",
        "vibration_energies_eV": "eV",
        "potential_energy_eV": "eV",
    },
    "magnetic-moments": {"species": "symbols", "moments_muB": "muB"},
    "elastic-moduli": {"stiffness_GPa": "GPa"},
    "hull": {"entries": "entries"},
    "eos-fit": {"volumes_angstrom3": "angstrom3", "energies_eV": "eV"},
    "intercalation-voltage": {"e_full": "eV", "e_empty": "eV", "e_metal": "eV"},
    "grain-boundary-energy": {
        "e_gb": "eV",
        "e_bulk_per_atom": "eV/atom",
        "area_A2": "angstrom2",
    },
    "error-metrics": {"predictions": "input-unit", "targets": "input-unit"},
    "spectrum-similarity": {
        "axis": "axis-unit",
        "first": "spectrum",
        "second": "spectrum",
    },
    "nmr-deconvolution": {
        "axis_ppm": "ppm",
        "mixture": "spectrum",
        "references": "spectrum",
    },
}
for _name in OPERATIONS:
    PORTS.setdefault(
        _name, {"text": "text"} if "text" in tools.SCHEMAS[_name]["properties"] else {}
    )


def selected(tool):
    if tool not in OPERATIONS:
        raise ValueError("Operation excluded or unknown in Anvil material catalog")
    enabled = os.environ.get("ANVIL_ENABLED_SKILL_TOOLS")
    if enabled is not None and tool not in enabled.split(","):
        raise ValueError("Operation is not enabled in this installation")
    return tool


def contract(tool):
    if tool not in OPERATIONS:
        raise ValueError("Unknown material operation")
    payload = dict(
        schema=SCHEMA,
        tool=tool,
        input_schema=tools.SCHEMAS[tool],
        ports=PORTS[tool],
        policy=CATALOG,
        publication="anvil.material-publication/v1",
    )
    # Contract bytes include output validation and mapping, separately from engine identity.
    payload["implementation"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(
        schema=SCHEMA,
        sha256=hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest(),
    )


def listing():
    return dict(
        **CATALOG,
        operations=[
            dict(
                tool=t,
                contract=contract(t),
                input_schema=tools.SCHEMAS[t],
                ports=PORTS[t],
                dependencies=tools.DEPENDENCIES[t],
                readiness=availability(t),
                outputs=dict(
                    scalars=SCALARS.get(t, []),
                    series=SERIES.get(t, []),
                    arrays=ARRAYS.get(t, []),
                    generated=t in GENERATED,
                    artifact_only=ARTIFACT_ONLY.get(t),
                ),
                enabled=os.environ.get("ANVIL_ENABLED_SKILL_TOOLS") is None
                or t in os.environ["ANVIL_ENABLED_SKILL_TOOLS"].split(","),
            )
            for t in OPERATIONS
        ],
    )


def validate_spec(spec):
    tool = selected(spec["tool"])
    expected = {"kind", "tool", "inputs", "implementation"}
    if "contract" in spec:
        expected |= {"contract", "bindings", "input_nodes"}
        if spec["contract"] != contract(tool):
            raise ValueError("Material provider contract changed")
        if not isinstance(spec["bindings"], list) or len(spec["bindings"]) > 2000:
            raise ValueError("Bounded explicit input bindings required")
    if set(spec) != expected or spec["kind"] != "skill-analysis":
        raise ValueError("Invalid skill specification")
    tools.validate(tool, spec["inputs"])


def validate_report(tool, report, inputs):
    """Additional v2 semantic checks; old provider report formats stay immutable."""
    tools.validate_report(tool, report)
    r = report["result"]
    for key in ("energy_basis", "response", "metric", "axis_unit"):
        if key in inputs and key in r and r[key] != inputs[key]:
            raise ValueError("Output metadata mismatch: " + key)
    if tool == "rdf" and (
        r["frame_count"] != len(inputs["frames"]) or len(r["g_r"]) != inputs["bins"]
    ):
        raise ValueError("RDF frame/bin count mismatch")
    if tool == "xrd":
        for families in r["hkls"]:
            if not isinstance(families, list) or not families:
                raise ValueError("Missing reflection family")
            for family in families:
                if (
                    set(family) != {"hkl", "multiplicity"}
                    or len(family["hkl"]) not in (3, 4)
                    or any(
                        not isinstance(i, int) or isinstance(i, bool)
                        for i in family["hkl"]
                    )
                    or not isinstance(family["multiplicity"], int)
                    or isinstance(family["multiplicity"], bool)
                    or family["multiplicity"] <= 0
                ):
                    raise ValueError("Invalid reflection family")

    if tool == "ideal-gas" and any(
        r[k] != inputs[k] for k in ("temperature_K", "pressure_Pa")
    ):
        raise ValueError("Thermodynamic conditions mismatch")
    if tool == "error-metrics" and (
        r["unit"] != inputs["unit"] or r["sample_count"] != len(inputs["targets"])
    ):
        raise ValueError("Metric units/count mismatch")
    if tool == "probability-density":
        import math
        import numpy as np

        lattice = np.asarray(inputs["frames"][0]["lattice"])
        expected = (
            np.ceil(np.linalg.norm(lattice, axis=1) / inputs["interval_angstrom"])
            .astype(int)
            .tolist()
        )
        if (
            not math.isclose(
                r["cell_volume_angstrom3"], np.linalg.det(lattice), rel_tol=1e-9
            )
            or r["shape"] != expected
        ):
            raise ValueError("Density cell/grid mismatch")
        shape = r["shape"]
        values = r["density_angstrom_minus3"]
        if len(values) != shape[0] or any(
            len(a) != shape[1] or any(len(b) != shape[2] for b in a) for a in values
        ):
            raise ValueError("Density shape mismatch")
        integral = (
            sum(v for a in values for b in a for v in b)
            * r["cell_volume_angstrom3"]
            / math.prod(shape)
        )
        if (
            not math.isclose(integral, 1, rel_tol=1e-6)
            or r["frame_count"] != len(inputs["frames"])
            or r["frame_interval_ps"] != inputs["frame_interval_ps"]
        ):
            raise ValueError("Density normalization/time mismatch")
    if tool == "magnetic-moments" and r["site_moments"] != inputs["moments_muB"]:
        raise ValueError("Site moment identity mismatch")


# Compound count is dimensionless; per-atom normalization remains in unit/code.
SCALARS = {
    "hull": [
        ("energy_above_hull_eV_atom", "eV/atom", "energy", ""),
        ("formation_energy_eV_atom", "eV/atom", "energy", ""),
    ],
    "eos-fit": [
        ("equilibrium_volume_angstrom3", "angstrom3", "volume", ""),
        ("equilibrium_energy_eV", "eV", "energy", "electronic"),
        ("bulk_modulus_GPa", "GPa", "pressure", ""),
    ],
    "ideal-gas": [
        ("enthalpy_eV", "eV", "energy", "enthalpy"),
        ("gibbs_energy_eV", "eV", "energy", "free_energy"),
        ("entropy_eV_K", "eV/K", "energy/temperature", ""),
    ],
    "orca-energy": [("final_energy_eV", "eV", "energy", "electronic")],
    "intercalation-voltage": [("voltage_V", "V", "electric_potential", "")],
    "grain-boundary-energy": [
        ("grain_boundary_energy_J_m2", "J/m2", "energy/area", "")
    ],
    "elastic-moduli": [
        (k, "GPa", "pressure", "")
        for k in ("bulk_modulus_GPa", "shear_modulus_GPa", "young_modulus_GPa")
    ]
    + [("poisson_ratio", "1", "dimensionless", "")],
    "error-metrics": [
        (k, "input-unit", "input-dimension", "") for k in ("mae", "rmse")
    ],
    "spectrum-similarity": [("score", "1", "dimensionless", "")],
    "magnetic-moments": [("total_magnetization", "muB", "magnetic_moment", "")],
    "symmetry": [("number", "1", "dimensionless", "")],
    "slab": [("surface_area_angstrom2", "angstrom2", "area", "")],
    "wannier-spreads": [
        (k, "angstrom2", "area", "")
        for k in ("omega_I_A2", "omega_D_A2", "omega_OD_A2", "omega_total_A2")
    ],
    "nmr-deconvolution": [
        ("wasserstein_distance", "ppm", "dimensionless", ""),
        ("noise", "1", "dimensionless", ""),
    ],
}
SERIES = {
    "rdf": [("r_angstrom", "angstrom", "length", "g_r", "1", "dimensionless")],
    "dielectric-chi": [
        ("energy_eV", "eV", "energy", k, "1", "dimensionless")
        for k in ("real", "imaginary")
    ],
    "edi-transport": [
        ("T_K", "K", "temperature", k, "cm2/V/s", "mobility")
        for k in ("mu_SERTA_xx", "mu_MRTA_xx", "mu_SERTA_yy", "mu_MRTA_yy")
    ],
}
ARRAYS = {
    "hull": [("decomposition", "1", "dimensionless")],
    "magnetic-moments": [("site_moments", "muB", "magnetic_moment")],
    "nmr-deconvolution": [("proportions", "1", "dimensionless")],
    "wannier-spreads": [
        ("wf_centres_xyz_A", "angstrom", "length"),
        ("wf_spreads_A2", "angstrom2", "area"),
    ],
    "edi-lifetimes": [
        ("E_eV", "eV", "energy"),
        ("inv_tau_SERTA_Ry", "Ry", "energy"),
        ("inv_tau_MRTA_Ry", "Ry", "energy"),
        ("tau_SERTA_fs", "fs", "time"),
        ("tau_MRTA_fs", "fs", "time"),
    ],
    "edi-matrix": [
        ("m2", "Ry2", "energy^2"),
        ("reM", "Ry", "energy"),
        ("imM", "Ry", "energy"),
        ("m2_loc", "Ry2", "energy^2"),
        ("m2_nl", "Ry2", "energy^2"),
    ],
    "epw-coupling-table": [
        ("enk", "eV", "energy"),
        ("enkq", "eV", "energy"),
        ("omega_meV", "meV", "energy"),
        ("g_meV", "meV", "energy"),
    ],
    "probability-density": [
        ("density_angstrom_minus3", "angstrom^-3", "inverse_volume")
    ],
}
GENERATED = {
    "supercell": None,
    "substitute": None,
    "slab": "structure",
    "symmetry": "conventional",
    "random-cubic-structure": "structure",
}
ARTIFACT_ONLY = {
    "uniform-kpoints": "Preprocessing reciprocal-fraction grid; not an observed property."
}


def density_axes(inputs, shape):
    # ProbabilityDensityAnalysis uses arange(0, 1, interval / lattice length),
    # which is not j / shape when interval does not divide the cell edge.
    import numpy as np

    lengths = np.linalg.norm(inputs["frames"][0]["lattice"], axis=1)
    return [
        dict(
            name="fractional_" + axis,
            kind="state",
            labels=[str(j * inputs["interval_angstrom"] / length) for j in range(n)],
            basis="Nearest grid-point occupancy on original lattice; equal weight cell volume/grid size, not physical Voronoi volume",
        )
        for axis, n, length in zip("abc", shape, lengths)
    ]


def availability(tool):
    from importlib.metadata import PackageNotFoundError

    missing = []
    for name in tools.DEPENDENCIES[tool]:
        try:
            tools.dependency_version(name)
        except PackageNotFoundError:
            missing.append(name)
    return (
        "Missing dependency metadata: " + ", ".join(missing)
        if missing
        else "Supported operation; dependency metadata present (not an execution certificate)"
    )


SUBJECT_PORTS = {
    "structure-match": "first",
    "nmr-deconvolution": "mixture",
    "spectrum-similarity": "first",
    "error-metrics": "predictions",
}
