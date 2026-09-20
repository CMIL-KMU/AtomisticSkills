"""Bounded analysis of explicit data using existing ASE and pymatgen APIs.

No calculator, checkpoint downloader, scheduler or arbitrary module dispatch.
Each result is provisional evidence, never a scientific acceptance decision.
"""

import hashlib
from importlib.metadata import version, PackageNotFoundError
import json
from pathlib import Path


def object_schema(**fields):
    return dict(
        type="object",
        properties=fields,
        required=list(fields),
        additionalProperties=False,
    )


def array(items, minimum=1, maximum=10000):
    return dict(type="array", items=items, minItems=minimum, maxItems=maximum)


NUMBER = dict(type="number")
POSITIVE = dict(type="number", exclusiveMinimum=0)
TEXT = dict(type="string", minLength=1, maxLength=256)
VECTOR = array(NUMBER, 3, 3)
STRUCTURE = object_schema(
    species=array(TEXT, maximum=2000),
    lattice=array(VECTOR, 3, 3),
    fractional=array(VECTOR, maximum=2000),
)
ENTRY = object_schema(identity=TEXT, formula=TEXT, energy_eV=NUMBER)
SCHEMAS = {
    "nmr-deconvolution": object_schema(
        axis_ppm=array(NUMBER, 2, 2000),
        mixture=array(dict(type="number", minimum=0), 2, 2000),
        references=array(array(dict(type="number", minimum=0), 2, 2000), 1, 20),
        protons=array(dict(type="integer", minimum=1), 1, 20),
        kappa_ppm=POSITIVE,
    ),
    "edi-transport": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "edi-lifetimes": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "edi-matrix": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "wannier-spreads": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "uniform-kpoints": object_schema(
        mesh=array(dict(type="integer", minimum=1, maximum=40), 3, 3)
    ),
    "epw-coupling-table": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "magnetic-moments": object_schema(
        species=array(TEXT, 1, 2000), moments_muB=array(NUMBER, 1, 2000)
    ),
    "dielectric-chi": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000),
        response=dict(enum=["ipa", "rpa"]),
    ),
    "probability-density": object_schema(
        frames=array(STRUCTURE, 1, 1000),
        species=TEXT,
        interval_angstrom=POSITIVE,
        frame_interval_ps=POSITIVE,
    ),
    "random-cubic-structure": object_schema(
        composition=TEXT,
        volume_scale=dict(type="number", minimum=0.5, maximum=4),
        seed=dict(type="integer", minimum=0, maximum=4294967295),
    ),
    "docking-box": object_schema(
        positions_angstrom=array(VECTOR, 1, 10000),
        padding_angstrom=dict(type="number", minimum=0, maximum=100),
        minimum_size_angstrom=POSITIVE,
    ),
    "redocking-rmsd": object_schema(
        reference_molblock=dict(type="string", minLength=1, maxLength=256000),
        probe_molblock=dict(type="string", minLength=1, maxLength=256000),
    ),
    "molecular-descriptors": object_schema(
        smiles=array(TEXT, 1, 250), include_sandp_tpsa=dict(type="boolean")
    ),
    "molecular-fingerprints": object_schema(
        smiles=array(TEXT, 1, 250),
        radius=dict(type="integer", minimum=1, maximum=6),
        fp_size=dict(enum=[512, 1024, 2048, 4096]),
        use_chirality=dict(type="boolean"),
        use_features=dict(type="boolean"),
    ),
    "orca-energy": object_schema(
        text=dict(type="string", minLength=1, maxLength=4000000)
    ),
    "error-metrics": object_schema(
        predictions=array(NUMBER),
        targets=array(NUMBER),
        unit=dict(enum=["eV/atom", "eV/angstrom", "eV/angstrom3"]),
        basis=TEXT,
    ),
    "spectrum-similarity": object_schema(
        axis=array(NUMBER, 2),
        first=array(NUMBER, 2),
        second=array(NUMBER, 2),
        axis_unit=dict(enum=["ppm", "cm^-1"]),
        metric=dict(enum=["l2", "cosine", "wasserstein"]),
    ),
    "intercalation-voltage": object_schema(
        e_full=NUMBER,
        e_empty=NUMBER,
        e_metal=NUMBER,
        n_metal=dict(type="integer", minimum=1),
        n_ions=dict(type="integer", minimum=1),
        metal_symbol=dict(enum=["Li", "Na", "K"]),
        energy_basis=TEXT,
    ),
    "grain-boundary-energy": object_schema(
        e_gb=NUMBER,
        n_atoms=dict(type="integer", minimum=1),
        e_bulk_per_atom=NUMBER,
        area_A2=POSITIVE,
        energy_basis=TEXT,
        interfaces=dict(const=2),
    ),
    "elastic-moduli": object_schema(stiffness_GPa=array(array(NUMBER, 6, 6), 6, 6)),
    "symmetry": object_schema(
        structure=STRUCTURE, symprec_angstrom=POSITIVE, angle_tolerance_degrees=POSITIVE
    ),
    "slab": object_schema(
        structure=STRUCTURE,
        miller=array(dict(type="integer", minimum=-6, maximum=6), 3, 3),
        thickness_angstrom=dict(type="number", exclusiveMinimum=0, maximum=100),
        vacuum_angstrom=dict(type="number", exclusiveMinimum=0, maximum=100),
        shift=dict(type="number", minimum=0, exclusiveMaximum=1),
    ),
    "xrd": object_schema(
        structure=STRUCTURE,
        wavelength_angstrom=POSITIVE,
        symprec_angstrom=dict(type="number", minimum=0),
        two_theta_degrees=array(NUMBER, 2, 2),
    ),
    "structure-match": object_schema(
        first=STRUCTURE,
        second=STRUCTURE,
        ltol=POSITIVE,
        stol=POSITIVE,
        angle_tol_degrees=POSITIVE,
    ),
    "supercell": object_schema(
        structure=STRUCTURE,
        repeats=array(dict(type="integer", minimum=1, maximum=20), 3, 3),
    ),
    "substitute": object_schema(
        structure=STRUCTURE,
        species_map=dict(
            type="object", minProperties=1, maxProperties=100, additionalProperties=TEXT
        ),
    ),
    "hull": object_schema(
        entries=array(ENTRY, 2, 2000), target=TEXT, energy_basis=TEXT
    ),
    "eos-fit": object_schema(
        volumes_angstrom3=array(POSITIVE, 5, 200),
        energies_eV=array(NUMBER, 5, 200),
        energy_basis=TEXT,
    ),
    "ideal-gas": object_schema(
        symbols=array(TEXT, 1, 1000),
        positions_angstrom=array(VECTOR, 1, 1000),
        vibration_energies_eV=array(POSITIVE, 0, 3000),
        potential_energy_eV=NUMBER,
        temperature_K=POSITIVE,
        pressure_Pa=POSITIVE,
        geometry=dict(enum=["monatomic", "linear", "nonlinear"]),
        symmetry_number=dict(type="integer", minimum=1),
        spin=dict(type="number", minimum=0, multipleOf=0.5),
    ),
    "rdf": object_schema(
        frames=array(STRUCTURE, 1, 1000),
        elements=array(dict(type="integer", minimum=1, maximum=118), 2, 2),
        rmax_angstrom=POSITIVE,
        bins=dict(type="integer", minimum=2, maximum=1000),
    ),
}
DEPENDENCIES = {name: ["numpy", "scipy", "ase", "pymatgen"] for name in SCHEMAS}
DEPENDENCIES.update(
    {
        "orca-energy": [],
        "intercalation-voltage": ["numpy"],
        "grain-boundary-energy": ["numpy"],
        "error-metrics": ["numpy"],
        "spectrum-similarity": ["numpy", "scipy"],
        "molecular-descriptors": ["numpy", "rdkit"],
        "molecular-fingerprints": ["numpy", "rdkit"],
    }
)
DEPENDENCIES.update({"docking-box": ["numpy"], "redocking-rmsd": ["numpy", "rdkit"]})
DEPENDENCIES.update(
    {
        "magnetic-moments": ["numpy"],
        "dielectric-chi": ["numpy"],
        "probability-density": [
            "numpy",
            "scipy",
            "pymatgen",
            "pymatgen-analysis-diffusion",
        ],
    }
)


DEPENDENCIES.update(
    {
        name: []
        for name in [
            "edi-transport",
            "edi-lifetimes",
            "edi-matrix",
            "wannier-spreads",
            "uniform-kpoints",
            "epw-coupling-table",
        ]
    }
)
DEPENDENCIES["nmr-deconvolution"] = ["numpy", "scipy"]


for _dependencies in DEPENDENCIES.values():
    if "pymatgen" in _dependencies:
        _dependencies.append("pymatgen-core")


def validate(tool, inputs):
    """Validate before importing optional engines; reject NaN and unknown keys."""
    from jsonschema import Draft202012Validator

    if tool not in SCHEMAS:
        raise ValueError("Unknown analysis tool")
    raw = json.dumps(inputs, allow_nan=False)
    if len(raw.encode()) > 8 * 1024 * 1024:
        raise ValueError("Analysis input exceeds the bounded inline contract")
    Draft202012Validator(SCHEMAS[tool]).validate(inputs)
    return inputs


def structure(value):
    from pymatgen.core import Structure
    import numpy as np

    if (
        len(value["species"]) != len(value["fractional"])
        or np.linalg.det(value["lattice"]) <= 0
    ):
        raise ValueError("Ordered sites and a positive-volume lattice required")
    result = Structure(value["lattice"], value["species"], value["fractional"])
    if not result.is_ordered:
        raise ValueError("Only ordered structures are supported")
    return result


def geometry(value):
    return dict(
        species=[str(s.specie) for s in value],
        lattice=value.lattice.matrix.tolist(),
        fractional=value.frac_coords.tolist(),
    )


def run(tool, inputs):
    """Call one audited library operation and return JSON data with explicit units."""
    p = validate(tool, inputs)
    if DEPENDENCIES[tool]:
        import numpy as np
    if tool == "nmr-deconvolution":
        from atomistic_analysis.nmr import deconvolve_spectra

        axis = np.asarray(p["axis_ppm"])
        rows = [np.asarray(p["mixture"]), *map(np.asarray, p["references"])]
        if (
            any(len(r) != len(axis) or r.sum() <= 0 for r in rows)
            or np.any(np.diff(np.round(axis, 6)) <= 0)
            or not np.allclose(np.diff(axis), np.diff(axis)[0])
            or len(p["protons"]) != len(p["references"])
        ):
            raise ValueError(
                "Paired positive-area spectra on one increasing uniform ppm grid and explicit proton counts required"
            )
        if np.linalg.matrix_rank(np.array(rows[1:])) != len(rows) - 1:
            raise ValueError("Reference spectra are linearly dependent")
        result = deconvolve_spectra(
            np.column_stack([axis, rows[0]]),
            [np.column_stack([axis, r]) for r in rows[1:]],
            p["protons"],
            kappa=p["kappa_ppm"],
        )
        if sum(result["proportions"]) <= 0:
            raise ValueError("No component signal recovered")
        result.update(
            scope="Supplied known-component proton-corrected mole fractions; baseline and quantitative NMR assumptions require review",
            distance_unit="ppm",
        )
    elif tool in ("edi-transport", "edi-lifetimes", "edi-matrix", "epw-coupling-table"):
        from tempfile import TemporaryDirectory
        from atomistic_analysis.edi_transport import parse_transport
        from atomistic_analysis.edi_lifetime import parse_inv_tau, summarize
        from atomistic_analysis.edi_matrix import parse_edmat
        from atomistic_analysis.epw import parse_rows

        rows = [
            line.split()
            for line in p["text"].splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        widths = {len(row) for row in rows}
        expected = {
            "edi-transport": {5},
            "edi-lifetimes": {7},
            "edi-matrix": {9, 13, 15},
            "epw-coupling-table": {7},
        }[tool]
        if len(widths) != 1 or not widths <= expected or len(rows) > 10000:
            raise ValueError(
                "One exact supported table layout with at most 10000 rows required"
            )
        with TemporaryDirectory() as directory:
            path = Path(directory) / "table.dat"
            path.write_text(p["text"])
            if tool == "edi-transport":
                data, metadata = parse_transport(str(path))
                if any(
                    r["T_K"] <= 0 or any(v < 0 for k, v in r.items() if k != "T_K")
                    for r in data
                ):
                    raise ValueError(
                        "Positive temperatures and nonnegative mobilities required"
                    )
                result = dict(rows=data, metadata=metadata, mobility_unit="cm2/V/s")
            elif tool == "edi-lifetimes":
                data = parse_inv_tau(str(path))
                if any(
                    r["ik"] < 1
                    or r["ibnd"] < 1
                    or any(
                        r[k] < 0
                        for k in (
                            "tau_SERTA_fs",
                            "tau_MRTA_fs",
                            "inv_tau_SERTA_Ry",
                            "inv_tau_MRTA_Ry",
                        )
                    )
                    for r in data
                ):
                    raise ValueError(
                        "Positive state indices and nonnegative lifetimes/rates required"
                    )
                result = dict(rows=data, summary=summarize(data))
            elif tool == "edi-matrix":
                data, layout = parse_edmat(str(path))
                result = dict(
                    rows=data,
                    columns=layout,
                    units=dict(m2="Ry2", reM="Ry", imM="Ry"),
                    scope="Printed EDI matrix values; caller retains provider normalization, no conversion or mobility inference",
                )
            else:
                data = parse_rows("\n".join(" ".join(r) for r in rows))
                if len(data) != len(rows):
                    raise ValueError("Unsupported coupling table numbers")
                result = dict(
                    rows=data,
                    scope="Selected printed coupling table only; no automatic CBM selection or top-mode gauge-invariance claim",
                )
        result["scientific_scope"] = (
            "Parsed supplied output only; no engine execution or convergence acceptance"
        )
    elif tool == "wannier-spreads":
        from atomistic_analysis.wannier import parse_wout
        import re

        result = parse_wout(p["text"])
        if (
            "Final State" not in p["text"]
            or not result["wf_spreads_A2"]
            or any(
                result[k] is None or result[k] < 0
                for k in ("omega_I_A2", "omega_D_A2", "omega_OD_A2", "omega_total_A2")
            )
            or re.search(r"WF centre and spread[^\n]*[0-9.][EeDd][+-]?\d", p["text"])
        ):
            raise ValueError("Complete supported final Wannier block required")
        result["scope"] = (
            "Printed centers and spreads; localization convergence and image folding require review"
        )
    elif tool == "uniform-kpoints":
        from atomistic_analysis.kpoints import uniform_kpoints

        if p["mesh"][0] * p["mesh"][1] * p["mesh"][2] > 10000:
            raise ValueError("K grid exceeds 10000 points")
        result = dict(
            fractional_reciprocal=[list(v) for v in uniform_kpoints(*p["mesh"])],
            scope="Uniform unshifted full grid; no symmetry reduction or electronic convergence claim",
        )
    elif tool == "magnetic-moments":
        from atomistic_analysis.magnetism import parse_magnetic_moments

        if len(p["species"]) != len(p["moments_muB"]):
            raise ValueError("Paired species and collinear moments required")
        result = parse_magnetic_moments(
            dict(
                data=dict(
                    magmom=p["moments_muB"],
                    structure=dict(
                        sites=[dict(species=[dict(element=s)]) for s in p["species"]]
                    ),
                )
            )
        )
        result.update(
            unit="muB",
            scope="Collinear projected site moments and threshold 0.1 muB heuristic; not an interstitial-inclusive total or magnetic ground state",
        )
    elif tool == "dielectric-chi":
        from atomistic_analysis.dielectric import parse_dielectric_section
        from tempfile import TemporaryDirectory

        header = {
            "ipa": "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)",
            "rpa": "INVERSE MACROSCOPIC DIELECTRIC TENSOR",
        }[p["response"]]
        if (
            p["text"].count(header) != 1
            or "HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR" in p["text"]
        ):
            raise ValueError("One explicit frequency-dependent CHI block required")
        with TemporaryDirectory() as directory:
            path = Path(directory) / "OUTCAR"
            path.write_text(p["text"])
            energy, real, imag = parse_dielectric_section(path, header)
        if len(energy) < 2 or np.any(np.diff(energy) <= 0):
            raise ValueError("Increasing photon energy samples required")
        result = dict(
            energy_eV=energy.tolist(),
            real=real.tolist(),
            imaginary=imag.tolist(),
            response=p["response"],
            scope="Scalar CHI block as printed; RPA inverse-header convention retained without inventing a tensor",
        )
    elif tool == "probability-density":
        from pymatgen.analysis.diffusion.aimd.pathway import ProbabilityDensityAnalysis

        frames = [structure(v) for v in p["frames"]]
        first = frames[0]
        if p["species"] not in [s.specie.symbol for s in first] or not np.allclose(
            first.lattice.angles, [90, 90, 90]
        ):
            raise ValueError(
                "Present mobile species and orthogonal reference cell required"
            )
        if any(
            f.species != first.species
            or not np.allclose(f.lattice.matrix, first.lattice.matrix)
            for f in frames
        ):
            raise ValueError("Constant lattice and site ordering required")
        if (
            np.prod(np.ceil(np.array(first.lattice.abc) / p["interval_angstrom"]))
            > 50000
        ):
            raise ValueError("Probability grid exceeds 50000 points")
        density = ProbabilityDensityAnalysis(
            first,
            np.array([f.frac_coords for f in frames]),
            interval=p["interval_angstrom"],
            species=[p["species"]],
        )
        result = dict(
            density_angstrom_minus3=density.Pr.tolist(),
            shape=list(density.lens),
            cell_volume_angstrom3=float(first.volume),
            frame_count=len(frames),
            frame_interval_ps=p["frame_interval_ps"],
            scope="Raw normalized occupancy in supplied coordinate frame; no smoothing, drift removal, or inferred conduction paths",
        )
    elif tool == "random-cubic-structure":
        from pymatgen.core import Composition
        from atomistic_analysis.random_structures import generate_random_structure

        composition = Composition(p["composition"])
        if not 1 <= composition.num_atoms <= 100 or any(
            v != int(v) or v <= 0 for v in composition.values()
        ):
            raise ValueError(
                "Positive integer composition with at most 100 sites required"
            )
        state = np.random.get_state()
        try:
            np.random.seed(p["seed"])
            generated = generate_random_structure(
                composition, 221, volume_scale=p["volume_scale"], max_attempts=50
            )
        finally:
            np.random.set_state(state)
        if generated is None:
            raise ValueError("No candidate passed provider distance screening")
        result = dict(
            structure=geometry(generated),
            scope="Seeded random coordinates in cubic metric; P1 candidate, no requested space-group symmetry or stability claim",
        )
    elif tool == "docking-box":
        from atomistic_analysis.coordinates import compute_box

        result = dict(
            **compute_box(
                np.asarray(p["positions_angstrom"]),
                p["padding_angstrom"],
                p["minimum_size_angstrom"],
            ),
            unit="angstrom",
            scope="Supplied selected coordinates only; no automatic binding site detection",
        )
    elif tool == "redocking-rmsd":
        from rdkit import Chem
        from atomistic_analysis.molecules import (
            check_molecule_identity,
            symmetry_corrected_rmsd,
        )

        pair = [
            Chem.MolFromMolBlock(p[key], removeHs=True)
            for key in ("reference_molblock", "probe_molblock")
        ]
        if any(
            m is None
            or not 1 <= m.GetNumAtoms() <= 200
            or m.GetNumConformers() != 1
            or not m.GetConformer().Is3D()
            for m in pair
        ):
            raise ValueError("Two bounded 3D molecular blocks required")
        if Chem.MolToSmiles(pair[0]) != Chem.MolToSmiles(pair[1]):
            raise ValueError(
                "Exact molecular identity, charge and stereochemistry required"
            )
        check_molecule_identity(*pair)
        result = dict(
            rmsd_angstrom=symmetry_corrected_rmsd(*pair),
            scope="Heavy-atom in-place RMSD; both poses must share the receptor coordinate frame",
        )
    elif tool in ("molecular-descriptors", "molecular-fingerprints"):
        from atomistic_analysis.molecules import (
            compute_descriptors,
            compute_fingerprints,
        )

        if tool == "molecular-descriptors":
            rows = [
                compute_descriptors(s, include_sandp_tpsa=p["include_sandp_tpsa"])
                for s in p["smiles"]
            ]
            if not all(row["valid"] for row in rows):
                raise ValueError("Invalid molecule or unavailable descriptor")
            result = dict(
                descriptors=rows,
                units=dict(molecular_weight="g/mol", tpsa="angstrom2"),
                scope="Molecular descriptors and heuristics; no ADMET measurement or clinical assessment",
            )
        else:
            result = compute_fingerprints(
                [(s, None) for s in p["smiles"]],
                **{k: v for k, v in p.items() if k != "smiles"},
            )
            if result["n_valid"] != len(p["smiles"]):
                raise ValueError("Invalid molecule or unavailable fingerprint")
    elif tool == "orca-energy":
        from atomistic_analysis.orca import parse_energy
        import re

        if "ORCA TERMINATED NORMALLY" not in p["text"] or re.search(
            r"FINAL SINGLE POINT ENERGY\s+[-\d.]+[EeDd]", p["text"]
        ):
            raise ValueError(
                "Normal ORCA termination and supported decimal energy format required"
            )
        result = parse_energy(p["text"])
        if "final_energy_eV" not in result:
            raise ValueError("No final ORCA energy found")
        result["scope"] = (
            "Parsed supplied output; normal termination alone does not establish SCF or scientific convergence"
        )
    elif tool == "error-metrics":
        from atomistic_analysis.metrics import evaluate_metrics

        if len(p["predictions"]) != len(p["targets"]):
            raise ValueError("Paired predictions and targets required")
        result = dict(
            **evaluate_metrics(p["predictions"], p["targets"]),
            unit=p["unit"],
            basis=p["basis"],
            sample_count=len(p["targets"]),
        )
    elif tool == "spectrum-similarity":
        from atomistic_analysis.spectra import (
            normalize,
            similarity_l2,
            similarity_cosine,
            similarity_wasserstein,
        )

        x, first, second = [np.asarray(p[key]) for key in ("axis", "first", "second")]
        if (
            len(x) != len(first)
            or len(x) != len(second)
            or np.any(np.diff(x) <= 0)
            or not np.allclose(np.diff(x), np.diff(x)[0])
            or np.ptp(first) == 0
            or np.ptp(second) == 0
        ):
            raise ValueError(
                "Nonflat spectra on one paired increasing uniform axis required"
            )
        methods = dict(
            l2=similarity_l2,
            cosine=similarity_cosine,
            wasserstein=similarity_wasserstein,
        )
        result = dict(
            score=methods[p["metric"]](normalize(first), normalize(second)),
            metric=p["metric"],
            axis_unit=p["axis_unit"],
            scope="Min-max normalized supplied spectra; no molecule identification or reference retrieval",
        )
    elif tool == "intercalation-voltage":
        from atomistic_analysis.scalars import calculate_voltage

        result = calculate_voltage(
            **{k: v for k, v in p.items() if k != "energy_basis"}
        )
        result.update(
            energy_basis=p["energy_basis"],
            scope="Monovalent transfer; common host stoichiometry and energy reference require review",
        )
    elif tool == "grain-boundary-energy":
        from atomistic_analysis.scalars import compute_gb_energy

        result = dict(
            grain_boundary_energy_J_m2=compute_gb_energy(
                **{
                    k: v
                    for k, v in p.items()
                    if k not in ("energy_basis", "interfaces")
                }
            ),
            energy_basis=p["energy_basis"],
            scope="Two equivalent interfaces; supplied area must be the actual cross-product area",
        )
    elif tool == "elastic-moduli":
        from pymatgen.analysis.elasticity.elastic import ElasticTensor

        matrix = np.array(p["stiffness_GPa"])
        if not np.allclose(matrix, matrix.T) or np.min(np.linalg.eigvalsh(matrix)) <= 0:
            raise ValueError("Symmetric positive definite stiffness required")
        tensor = ElasticTensor.from_voigt(matrix)
        result = dict(
            bulk_modulus_GPa=float(tensor.k_vrh),
            shear_modulus_GPa=float(tensor.g_vrh),
            young_modulus_GPa=float(tensor.y_mod / 1e9),
            poisson_ratio=float(tensor.homogeneous_poisson),
            scope="Voigt-Reuss-Hill aggregate of supplied stiffness; no strain calculations",
        )
    elif tool == "symmetry":
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        analyzer = SpacegroupAnalyzer(
            structure(p["structure"]),
            symprec=p["symprec_angstrom"],
            angle_tolerance=p["angle_tolerance_degrees"],
        )
        result = dict(
            number=analyzer.get_space_group_number(),
            symbol=analyzer.get_space_group_symbol(),
            conventional=geometry(analyzer.get_conventional_standard_structure()),
        )
    elif tool == "slab":
        from pymatgen.core.surface import SlabGenerator

        if not any(p["miller"]):
            raise ValueError("Nonzero Miller index required")
        slab = SlabGenerator(
            structure(p["structure"]),
            p["miller"],
            p["thickness_angstrom"],
            p["vacuum_angstrom"],
            center_slab=True,
        ).get_slab(shift=p["shift"])
        if len(slab) > 2000:
            raise ValueError("Slab exceeds 2000 sites")
        result = dict(
            structure=geometry(slab),
            surface_area_angstrom2=float(slab.surface_area),
            scope="Unrelaxed termination; polarity, stoichiometry and convergence require review",
        )
    elif tool == "xrd":
        from pymatgen.analysis.diffraction.xrd import XRDCalculator

        low, high = p["two_theta_degrees"]
        if not 0 <= low < high <= 180:
            raise ValueError("Two-theta range must increase within [0, 180] degrees")
        pattern = XRDCalculator(
            wavelength=p["wavelength_angstrom"], symprec=p["symprec_angstrom"]
        ).get_pattern(structure(p["structure"]), two_theta_range=(low, high))
        result = dict(
            two_theta_degrees=pattern.x.tolist(),
            intensity_relative=pattern.y.tolist(),
            d_angstrom=list(pattern.d_hkls),
            hkls=pattern.hkls,
        )
    elif tool == "structure-match":
        from pymatgen.analysis.structure_matcher import StructureMatcher

        match = StructureMatcher(
            ltol=p["ltol"], stol=p["stol"], angle_tol=p["angle_tol_degrees"]
        )
        result = dict(
            matches=bool(match.fit(structure(p["first"]), structure(p["second"]))),
            scope="Provided structures only; no database or literature novelty claim",
        )
    elif tool in ("supercell", "substitute"):
        value = structure(p["structure"])
        if tool == "supercell":
            if len(value) * np.prod(p["repeats"]) > 2000:
                raise ValueError("Supercell exceeds 2000 sites")
            value.make_supercell(p["repeats"])
        else:
            if not set(p["species_map"]) <= set(p["structure"]["species"]):
                raise ValueError("Substitution species absent from structure")
            value.replace_species(p["species_map"])
        result = geometry(value)
    elif tool == "hull":
        from pymatgen.analysis.phase_diagram import PhaseDiagram
        from pymatgen.entries.computed_entries import ComputedEntry

        if len({r["identity"] for r in p["entries"]}) != len(p["entries"]):
            raise ValueError("Entry identities must be unique")
        entries = [
            ComputedEntry(r["formula"], r["energy_eV"], entry_id=r["identity"])
            for r in p["entries"]
        ]
        matches = [e for e in entries if e.entry_id == p["target"]]
        if len(matches) != 1:
            raise ValueError("Explicit target entry identity required")
        pd = PhaseDiagram(entries)
        decomposition, energy = pd.get_decomp_and_e_above_hull(matches[0])
        result = dict(
            energy_above_hull_eV_atom=float(energy),
            formation_energy_eV_atom=float(pd.get_form_energy_per_atom(matches[0])),
            decomposition=[
                dict(identity=e.entry_id, fraction=float(v))
                for e, v in decomposition.items()
            ],
            energy_basis=p["energy_basis"],
            scope="Provided phases only; completeness and common energy reference require review",
        )
    elif tool == "eos-fit":
        from ase.eos import EquationOfState
        from ase.units import GPa

        v, e = p["volumes_angstrom3"], p["energies_eV"]
        if len(v) != len(e) or len(set(v)) != len(v):
            raise ValueError("Paired distinct volume/energy samples required")
        volume, energy, bulk = EquationOfState(v, e, eos="birchmurnaghan").fit()
        if not min(v) < volume < max(v) or bulk <= 0:
            raise ValueError("Fit has no physical minimum inside the supplied scan")
        result = dict(
            equilibrium_volume_angstrom3=float(volume),
            equilibrium_energy_eV=float(energy),
            bulk_modulus_GPa=float(bulk / GPa),
            energy_basis=p["energy_basis"],
        )
    elif tool == "ideal-gas":
        from ase import Atoms
        from ase.thermochemistry import IdealGasThermo

        if len(p["symbols"]) != len(p["positions_angstrom"]):
            raise ValueError("Paired symbols/positions required")
        atoms = Atoms(p["symbols"], positions=p["positions_angstrom"])
        rank = np.linalg.matrix_rank(atoms.positions - atoms.positions[0])
        if (p["geometry"] == "linear" and rank != 1) or (
            p["geometry"] == "nonlinear" and rank < 2
        ):
            raise ValueError("Coordinates do not match the stated molecular geometry")
        expected = {
            "monatomic": 0,
            "linear": 3 * len(atoms) - 5,
            "nonlinear": 3 * len(atoms) - 6,
        }[p["geometry"]]
        if (
            expected < 0
            or len(p["vibration_energies_eV"]) != expected
            or (p["geometry"] == "monatomic" and len(atoms) != 1)
        ):
            raise ValueError(
                "Supply only the complete physical vibrational modes for the stated geometry"
            )
        thermo = IdealGasThermo(
            p["vibration_energies_eV"],
            p["geometry"],
            potentialenergy=p["potential_energy_eV"],
            atoms=atoms,
            symmetrynumber=p["symmetry_number"],
            spin=p["spin"],
            ignore_imag_modes=False,
        )
        t, pressure = p["temperature_K"], p["pressure_Pa"]
        result = dict(
            enthalpy_eV=thermo.get_enthalpy(t, verbose=False),
            entropy_eV_K=thermo.get_entropy(t, pressure, verbose=False),
            gibbs_energy_eV=thermo.get_gibbs_energy(t, pressure, verbose=False),
            temperature_K=t,
            pressure_Pa=pressure,
            approximation="Ideal gas, rigid rotor, harmonic oscillator; spin is S, not multiplicity",
        )
    else:
        from ase.geometry.rdf import get_rdf
        from pymatgen.io.ase import AseAtomsAdaptor

        frames = [AseAtomsAdaptor.get_atoms(structure(v)) for v in p["frames"]]
        for frame in frames:
            if not set(p["elements"]) <= set(frame.numbers):
                raise ValueError("RDF species absent from a frame")
        rows = [
            get_rdf(a, p["rmax_angstrom"], p["bins"], elements=p["elements"])
            for a in frames
        ]
        result = dict(
            r_angstrom=rows[0][1].tolist(),
            g_r=np.mean([r[0] for r in rows], axis=0).tolist(),
            frame_count=len(frames),
            normalization="ASE partial RDF; no coordination or diffusion convergence claim",
        )
    json.dumps(result, allow_nan=False)
    return dict(status="provisional", tool=tool, result=result)


def validate_report(tool, report):
    """Collection validates typed evidence without rerunning numerical engines."""
    from jsonschema import Draft202012Validator

    nonnegative = dict(type="number", minimum=0)
    scalar_fields = {
        "intercalation-voltage": ["voltage_V"],
        "grain-boundary-energy": ["grain_boundary_energy_J_m2"],
        "elastic-moduli": [
            "bulk_modulus_GPa",
            "shear_modulus_GPa",
            "young_modulus_GPa",
            "poisson_ratio",
        ],
        "eos-fit": [
            "equilibrium_volume_angstrom3",
            "equilibrium_energy_eV",
            "bulk_modulus_GPa",
        ],
        "hull": ["energy_above_hull_eV_atom", "formation_energy_eV_atom"],
        "ideal-gas": ["enthalpy_eV", "entropy_eV_K", "gibbs_energy_eV"],
        "orca-energy": ["final_energy_hartree", "final_energy_eV"],
        "error-metrics": ["mae", "rmse"],
    }
    required = {
        name: {key: NUMBER for key in fields} for name, fields in scalar_fields.items()
    }
    required.update(
        {
            "nmr-deconvolution": dict(
                proportions=array(nonnegative, 1, 20),
                wasserstein_distance=NUMBER,
                noise=nonnegative,
            ),
            **{
                name: dict(rows=array(dict(type="object")))
                for name in (
                    "edi-transport",
                    "edi-lifetimes",
                    "edi-matrix",
                    "epw-coupling-table",
                )
            },
            "wannier-spreads": dict(
                omega_total_A2=nonnegative,
                wf_centres_xyz_A=array(VECTOR),
                wf_spreads_A2=array(nonnegative),
            ),
            "uniform-kpoints": dict(fractional_reciprocal=array(VECTOR)),
            "magnetic-moments": dict(
                total_magnetization=NUMBER,
                site_moments=array(NUMBER),
                magnetic_ordering=TEXT,
            ),
            "dielectric-chi": dict(
                energy_eV=array(NUMBER, 2),
                real=array(NUMBER, 2),
                imaginary=array(NUMBER, 2),
            ),
            "probability-density": dict(
                density_angstrom_minus3=array(array(array(nonnegative))),
                shape=array(dict(type="integer", minimum=1), 3, 3),
                cell_volume_angstrom3=POSITIVE,
            ),
            "random-cubic-structure": dict(structure=STRUCTURE),
            "docking-box": {
                key: NUMBER
                for key in (
                    "center_x",
                    "center_y",
                    "center_z",
                    "size_x",
                    "size_y",
                    "size_z",
                )
            },
            "redocking-rmsd": dict(rmsd_angstrom=nonnegative),
            "xrd": dict(
                two_theta_degrees=array(NUMBER, 0),
                intensity_relative=array(nonnegative, 0),
                d_angstrom=array(POSITIVE, 0),
                hkls=dict(type="array"),
            ),
            "supercell": STRUCTURE["properties"],
            "substitute": STRUCTURE["properties"],
            "structure-match": dict(matches=dict(type="boolean")),
            "rdf": dict(
                r_angstrom=array(POSITIVE),
                g_r=array(nonnegative),
                frame_count=dict(type="integer", minimum=1),
            ),
            "symmetry": dict(
                number=dict(type="integer", minimum=1, maximum=230),
                symbol=TEXT,
                conventional=STRUCTURE,
            ),
            "slab": dict(structure=STRUCTURE, surface_area_angstrom2=POSITIVE),
            "spectrum-similarity": dict(
                score=dict(type="number", minimum=0, maximum=1)
            ),
            "molecular-descriptors": dict(
                descriptors=array(dict(type="object"), 1, 250)
            ),
            "molecular-fingerprints": dict(
                n_compounds=dict(type="integer", minimum=1),
                n_valid=dict(type="integer", minimum=1),
                compounds=array(dict(type="object"), 1, 250),
                similarity_matrix=array(array(nonnegative), 1, 250),
            ),
        }
    )
    fields = required[tool]
    schema = object_schema(
        status=dict(const="provisional"),
        tool=dict(const=tool),
        result=dict(type="object", properties=fields, required=list(fields)),
    )
    json.dumps(report, allow_nan=False)
    Draft202012Validator(schema).validate(report)
    result = report["result"]
    paired = {
        "xrd": ["two_theta_degrees", "intensity_relative", "d_angstrom", "hkls"],
        "rdf": ["r_angstrom", "g_r"],
        "supercell": ["species", "fractional"],
        "substitute": ["species", "fractional"],
        "dielectric-chi": ["energy_eV", "real", "imaginary"],
    }
    if tool in paired and len({len(result[k]) for k in paired[tool]}) != 1:
        raise ValueError("Mismatched scientific output arrays")


def dependency_version(name):
    """New pymatgen releases split core; older releases bundle it in pymatgen."""
    try:
        return version(name)
    except PackageNotFoundError:
        if name == "pymatgen-core":
            return "bundled-in-pymatgen:" + version("pymatgen")
        raise


def identity(tool):
    """Science bytes and installed numerical dependencies are distinct evidence."""
    sources = b"".join(
        Path(__file__).with_name(name).read_bytes()
        for name in (
            "tools.py",
            "scalars.py",
            "orca.py",
            "spectra.py",
            "metrics.py",
            "molecules.py",
            "coordinates.py",
            "nmr.py",
            "edi_transport.py",
            "edi_lifetime.py",
            "edi_matrix.py",
            "wannier.py",
            "kpoints.py",
            "epw.py",
            "magnetism.py",
            "dielectric.py",
            "random_structures.py",
        )
    )
    return dict(
        source_sha256=hashlib.sha256(sources).hexdigest(),
        versions={
            n: dependency_version(n)
            for n in ["atomistic-analysis", "jsonschema", *DEPENDENCIES[tool]]
        },
    )
