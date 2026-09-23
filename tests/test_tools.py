"""Synthetic references for bounded library integrations; no model or network."""

import math
import pytest
from src.utils.analysis import tools


@pytest.fixture
def cube():
    return dict(
        species=["Na", "Cl"],
        lattice=[[6, 0, 0], [0, 6, 0], [0, 0, 6]],
        fractional=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )


def test_supercell_substitution_and_comparison(cube):
    expanded = tools.run("supercell", dict(structure=cube, repeats=[2, 1, 1]))["result"]
    assert len(expanded["species"]) == 4 and expanded["lattice"][0][0] == 12
    changed = tools.run("substitute", dict(structure=cube, species_map={"Na": "K"}))[
        "result"
    ]
    assert (
        changed["species"] == ["K", "Cl"]
        and changed["fractional"] == cube["fractional"]
    )
    p = dict(first=cube, second=expanded, ltol=0.2, stol=0.3, angle_tol_degrees=5)
    assert tools.run("structure-match", p)["result"]["matches"]
    p["second"] = changed
    assert not tools.run("structure-match", p)["result"]["matches"]


def test_xrd_simple_cubic_bragg_reference():
    cube = dict(
        species=["Si"],
        lattice=[[4, 0, 0], [0, 4, 0], [0, 0, 4]],
        fractional=[[0, 0, 0]],
    )
    r = tools.run(
        "xrd",
        dict(
            structure=cube,
            wavelength_angstrom=1,
            symprec_angstrom=0,
            two_theta_degrees=[0, 60],
        ),
    )["result"]
    assert r["two_theta_degrees"][0] == pytest.approx(
        2 * math.degrees(math.asin(1 / 8))
    )
    assert r["d_angstrom"][0] == pytest.approx(4)
    assert len(r["hkls"]) == len(r["two_theta_degrees"])


def test_hull_identity_not_formula_and_energy_reference():
    entries = [
        dict(identity=i, formula=f, energy_eV=e)
        for i, f, e in [
            ("li", "Li", 0),
            ("cl", "Cl", 0),
            ("low", "LiCl", -2),
            ("high", "LiCl", -1),
        ]
    ]
    r = tools.run(
        "hull",
        dict(entries=entries, target="high", energy_basis="synthetic common reference"),
    )["result"]
    assert r["energy_above_hull_eV_atom"] == pytest.approx(0.5)
    assert r["decomposition"] == [dict(identity="low", fraction=1)]
    with pytest.raises(ValueError):
        tools.run(
            "hull", dict(entries=entries, target="LiCl", energy_basis="synthetic")
        )


def test_ideal_gas_monatomic_reference():
    from ase.units import kB

    p = dict(
        symbols=["He"],
        positions_angstrom=[[0, 0, 0]],
        vibration_energies_eV=[],
        potential_energy_eV=-1,
        temperature_K=300,
        pressure_Pa=101325,
        geometry="monatomic",
        symmetry_number=1,
        spin=0,
    )
    r = tools.run("ideal-gas", p)["result"]
    assert r["enthalpy_eV"] == pytest.approx(-1 + 2.5 * kB * 300)
    assert r["gibbs_energy_eV"] == pytest.approx(
        r["enthalpy_eV"] - 300 * r["entropy_eV_K"]
    )
    with pytest.raises(ValueError):
        tools.run("ideal-gas", dict(p, vibration_energies_eV=[0.1]))
    with pytest.raises(ValueError, match="geometry"):
        tools.run(
            "ideal-gas",
            dict(
                p,
                symbols=["H", "H"],
                positions_angstrom=[[0, 0, 0], [0, 0, 0]],
                geometry="linear",
                vibration_energies_eV=[0.1],
            ),
        )


def test_rdf_pair_peak_and_small_cell_rejection(cube):
    cube["fractional"][1] = [1 / 6, 0, 0]
    r = tools.run(
        "rdf", dict(frames=[cube, cube], elements=[11, 17], rmax_angstrom=2, bins=20)
    )["result"]
    assert r["frame_count"] == 2 and sum(v > 0 for v in r["g_r"]) == 1
    assert r["r_angstrom"][max(range(20), key=lambda i: r["g_r"][i])] == pytest.approx(
        0.95
    )
    from ase.geometry.rdf import CellTooSmall

    with pytest.raises(CellTooSmall):
        tools.run(
            "rdf", dict(frames=[cube], elements=[11, 17], rmax_angstrom=8, bins=20)
        )


@pytest.mark.parametrize(
    "tool,inputs",
    [
        ("os.system", {}),
        ("xrd", {"unknown": 1}),
        ("hull", {"entries": []}),
        ("eos-fit", {"volumes_angstrom3": [float("nan")]}),
    ],
)
def test_invalid_inputs_fail_closed(tool, inputs):
    with pytest.raises(Exception):
        tools.run(tool, inputs)


def test_scalar_reference_and_tensor_units():
    r = tools.run(
        "intercalation-voltage",
        dict(
            e_full=-12,
            e_empty=-8,
            e_metal=-2,
            n_metal=2,
            n_ions=2,
            metal_symbol="Li",
            energy_basis="same method",
        ),
    )["result"]
    assert r["voltage_V"] == pytest.approx(1)
    r = tools.run(
        "grain-boundary-energy",
        dict(
            e_gb=-18,
            n_atoms=10,
            e_bulk_per_atom=-2,
            area_A2=10,
            interfaces=2,
            energy_basis="same method",
        ),
    )["result"]
    assert r["grain_boundary_energy_J_m2"] == pytest.approx(1.60218)
    c = [
        [
            120 if i == j and i < 3 else 40 if i < 3 and j < 3 else 40 if i == j else 0
            for j in range(6)
        ]
        for i in range(6)
    ]
    r = tools.run("elastic-moduli", dict(stiffness_GPa=c))["result"]
    assert r["bulk_modulus_GPa"] == pytest.approx(200 / 3)
    assert r["shear_modulus_GPa"] == pytest.approx(40)
    assert r["young_modulus_GPa"] == pytest.approx(100)
    c[0][0] = -1
    with pytest.raises(ValueError):
        tools.run("elastic-moduli", dict(stiffness_GPa=c))


def test_eos_independent_analytic_samples():
    # B'=4 Birch-Murnaghan analytic reference, V0=16, E0=-5, B0=.5 eV/A^3.
    volumes = [13, 14, 15, 16, 17, 18, 19]
    energies = [
        -5 + 9 * 16 * 0.5 / 16 * ((16 / v) ** (2 / 3) - 1) ** 2 * 2 for v in volumes
    ]
    r = tools.run(
        "eos-fit",
        dict(
            volumes_angstrom3=volumes,
            energies_eV=energies,
            energy_basis="synthetic cell",
        ),
    )["result"]
    assert r["equilibrium_volume_angstrom3"] == pytest.approx(16, rel=1e-6)
    assert r["equilibrium_energy_eV"] == pytest.approx(-5, rel=1e-6)
    assert r["bulk_modulus_GPa"] == pytest.approx(80.1088, rel=1e-5)


def test_symmetry_and_slab_geometry(cube):
    r = tools.run(
        "symmetry",
        dict(structure=cube, symprec_angstrom=0.01, angle_tolerance_degrees=5),
    )["result"]
    assert r["number"] == 221
    r = tools.run(
        "slab",
        dict(
            structure=cube,
            miller=[0, 0, 1],
            thickness_angstrom=8,
            vacuum_angstrom=12,
            shift=0,
        ),
    )["result"]
    assert r["surface_area_angstrom2"] == pytest.approx(36)
    assert r["structure"]["lattice"][2][2] >= 20


def test_optional_dependency_failure_is_not_success(monkeypatch, cube):
    import sys

    monkeypatch.setitem(sys.modules, "pymatgen.symmetry.analyzer", None)
    with pytest.raises(ImportError):
        tools.run(
            "symmetry",
            dict(structure=cube, symprec_angstrom=0.01, angle_tolerance_degrees=5),
        )


def test_existing_parsers_and_metrics():
    r = tools.run(
        "orca-energy",
        dict(text="FINAL SINGLE POINT ENERGY -1.0\nORCA TERMINATED NORMALLY"),
    )["result"]
    assert r["final_energy_eV"] == pytest.approx(-27.211386245988)
    for text in [
        "FINAL SINGLE POINT ENERGY -1.0",
        "FINAL SINGLE POINT ENERGY -1.0E-3\nORCA TERMINATED NORMALLY",
    ]:
        with pytest.raises(ValueError):
            tools.run("orca-energy", dict(text=text))
    r = tools.run(
        "error-metrics",
        dict(predictions=[1, 3], targets=[0, 1], unit="eV/atom", basis="paired labels"),
    )["result"]
    assert r["mae"] == pytest.approx(1.5) and r["rmse"] == pytest.approx(math.sqrt(2.5))
    for metric in ["l2", "cosine", "wasserstein"]:
        p = dict(
            axis=[0, 1, 2],
            first=[0, 2, 0],
            second=[0, 3, 0],
            axis_unit="ppm",
            metric=metric,
        )
        report = tools.run("spectrum-similarity", p)
        assert report["result"]["score"] == pytest.approx(1)
        tools.validate_report("spectrum-similarity", report)
    with pytest.raises(Exception):
        tools.validate_report(
            "hull", dict(status="provisional", tool="hull", result={})
        )


def test_molecule_utilities_without_partial_success():
    pytest.importorskip("rdkit")
    report = tools.run(
        "molecular-descriptors", dict(smiles=["CCO"], include_sandp_tpsa=False)
    )
    row = report["result"]["descriptors"][0]
    assert (
        row["molecular_weight"] == pytest.approx(46.07)
        and row["hbd"] == 1
        and row["hba"] == 1
    )
    p = dict(
        smiles=["CCO", "CCO"],
        radius=2,
        fp_size=512,
        use_chirality=True,
        use_features=False,
    )
    report = tools.run("molecular-fingerprints", p)
    assert report["result"]["similarity_matrix"] == [[1, 1], [1, 1]]
    tools.validate_report("molecular-fingerprints", report)
    with pytest.raises(ValueError):
        tools.run("molecular-fingerprints", dict(p, smiles=["CCO", "invalid"]))


def test_box_reference():
    r = tools.run(
        "docking-box",
        dict(
            positions_angstrom=[[0, 2, 4], [4, 4, 6]],
            padding_angstrom=2,
            minimum_size_angstrom=7,
        ),
    )["result"]
    assert [
        r[k] for k in ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    ] == [2, 3, 5, 8, 7, 7]


def test_redocking_does_not_align_away_translation():
    Chem = pytest.importorskip("rdkit.Chem")
    mol = Chem.MolFromSmiles("CO")
    conf = Chem.Conformer(2)
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0, 0, 0))
    conf.SetAtomPosition(1, (1.4, 0, 0))
    mol.AddConformer(conf)
    reference = Chem.MolToMolBlock(mol)
    for i in range(2):
        p = mol.GetConformer().GetAtomPosition(i)
        mol.GetConformer().SetAtomPosition(i, (p.x + 3, p.y, p.z))
    r = tools.run(
        "redocking-rmsd",
        dict(reference_molblock=reference, probe_molblock=Chem.MolToMolBlock(mol)),
    )
    assert r["result"]["rmsd_angstrom"] == pytest.approx(3)
    tools.validate_report("redocking-rmsd", r)


def test_magnetic_sign_and_frequency_dependent_dielectric():
    for moments in [[2, 2], [-2, -2]]:
        r = tools.run(
            "magnetic-moments", dict(species=["Fe", "Fe"], moments_muB=moments)
        )
        assert r["result"]["magnetic_ordering"] == "ferromagnetic"
        tools.validate_report("magnetic-moments", r)
    text = "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)\n0.0 2.0 0.0\n1.0 3.0 0.5\nend"
    r = tools.run("dielectric-chi", dict(text=text, response="ipa"))
    assert r["result"]["energy_eV"] == [0, 1] and r["result"]["imaginary"] == [0, 0.5]
    tools.validate_report("dielectric-chi", r)
    with pytest.raises(ValueError):
        tools.run(
            "dielectric-chi",
            dict(
                text=text + "\nHEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR",
                response="ipa",
            ),
        )


def test_density_normalization_and_fixed_cell(cube):
    pytest.importorskip("pymatgen.analysis.diffusion")
    import numpy as np

    p = dict(
        frames=[cube, cube], species="Na", interval_angstrom=2, frame_interval_ps=0.1
    )
    r = tools.run("probability-density", p)
    grid = np.array(r["result"]["density_angstrom_minus3"])
    assert grid.shape == (3, 3, 3) and np.count_nonzero(grid) == 1
    assert grid.sum() * 216 / grid.size == pytest.approx(1)
    tools.validate_report("probability-density", r)
    with pytest.raises(ValueError):
        tools.run("probability-density", dict(p, interval_angstrom=0.0001))


def test_seeded_cubic_candidate_is_not_claimed_symmetric():
    p = dict(composition="NaCl", volume_scale=2, seed=42)
    a = tools.run("random-cubic-structure", p)
    assert a == tools.run("random-cubic-structure", p)
    assert a["result"]["structure"]["species"] == ["Na", "Cl"]
    assert "P1 candidate" in a["result"]["scope"]
    with pytest.raises(ValueError):
        tools.run("random-cubic-structure", dict(p, composition="Na0.5Cl"))


def test_nmr_known_mixture_and_degenerate_reference_rejection():
    p = dict(
        axis_ppm=[0, 1, 2],
        mixture=[1, 0, 3],
        references=[[1, 0, 0], [0, 0, 1]],
        protons=[1, 3],
        kappa_ppm=2,
    )
    r = tools.run("nmr-deconvolution", p)
    assert r["result"]["proportions"] == pytest.approx([0.5, 0.5])
    assert r["result"]["wasserstein_distance"] == pytest.approx(0)
    tools.validate_report("nmr-deconvolution", r)
    with pytest.raises(ValueError, match="dependent"):
        tools.run("nmr-deconvolution", dict(p, references=[[1, 0, 0], [1, 0, 0]]))


def test_electronic_output_parsers_and_kgrid():
    r = tools.run(
        "edi-transport", dict(text="# Grid 2 2 1\n300 10 11 12 13\n600 5 6 7 8")
    )
    assert r["result"]["rows"][1]["mu_MRTA_yy"] == 8
    tools.validate_report("edi-transport", r)
    with pytest.raises(ValueError):
        tools.run("edi-transport", dict(text="300 1 2 3"))
    r = tools.run(
        "edi-lifetimes", dict(text="1 2 -3 0.1 0.2 10 20\n2 2 -2 0.2 0.3 20 40")
    )
    assert r["result"]["summary"]["mean_tau_SERTA_fs"] == 15
    r = tools.run("edi-matrix", dict(text="1 0 0 0 2 3 25 3 4"))
    assert r["result"]["rows"][0]["m2"] == 25 and r["result"]["columns"] == 9
    r = tools.run("epw-coupling-table", dict(text="2 2 1 -3 -3 12 4"))
    assert r["result"]["rows"][0]["g_meV"] == 4
    r = tools.run("uniform-kpoints", dict(mesh=[2, 1, 1]))
    assert r["result"]["fractional_reciprocal"] == [[0, 0, 0], [0.5, 0, 0]]
    tools.validate_report("uniform-kpoints", r)
    with pytest.raises(ValueError):
        tools.run("uniform-kpoints", dict(mesh=[40, 40, 40]))
    text = "Final State\nWF centre and spread 1 (0.0, 1.0, 2.0) 3.0\nOmega I = 1.0\nOmega D = 1.0\nOmega OD = 1.0\nOmega Total = 3.0"
    r = tools.run("wannier-spreads", dict(text=text))
    assert r["result"]["omega_total_A2"] == 3 and r["result"]["wf_centres_z_max"] == 2
    tools.validate_report("wannier-spreads", r)
    with pytest.raises(ValueError):
        tools.run("wannier-spreads", dict(text="Final State\nmissing"))
