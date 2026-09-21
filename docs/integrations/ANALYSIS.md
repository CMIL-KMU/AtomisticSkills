# Supplied-data analysis qualification

Active Anvil material policy and v2 DB contracts: [MATERIAL_DB.md](MATERIAL_DB.md). The provider functions below remain available upstream; only the explicit 27-operation material catalog is admitted by Anvil.

The 130-skill audit is in [MATRIX.md](MATRIX.md) and [skill-inventory.json](skill-inventory.json).
Each row retains its full source instructions, CLI/function contracts, units,
dependencies and exact remaining scope. These tutorials are not 130 executable
recipes. The small installed `atomistic_analysis.tools` registry exposes 31 bounded
operations; MLIP simulation remains Peregrine's responsibility and Anvil owns campaigns.
Native VASP DFT requests, including VASP MD/NEB, belong to the external engine
recipe lifecycle; they are not prohibited by the separate MLIP engine policy.

HTVS VASP integration uses `anvil-recipes`' `vasp.htvs.protocol` entrypoint and
`htvs/<legacy directory>/v1` scientific profiles. Anvil owns registered input
geometries, requester/attempt state and typed Calc/Result publication; mkite owns
native stage execution/parsing and `anvil-machine` owns site/Slurm capabilities.
Reuse the existing analysis v2 Geom/Result bindings for postprocessing, retaining
ResultTarget/ResultSource/ResultLink lineage. This package adds no VASP runner or
scheduler. The complete 25-profile migration matrix and explicit blocked families
are in `anvil-recipes/mkite_vaspsol/protocols/htvs_v1/PROTOCOL.md`; operator steps
are in `anvil/packages/anvil-db/anvil_mcp/db/VASP_WORKFLOW.md` in the sibling
workspace. Synthetic lifecycle tests do not establish physical convergence.

[operations.json](operations.json) contains exact allowlisted input schemas and
per-operation dependency identities. Inputs/outputs are finite JSON; unknown keys,
unknown operations, malformed arrays and oversized inline inputs are rejected.
All outputs remain provisional. No licensed binary, checkpoint, credential, remote
database or automatically chosen model is silently downloaded or substituted.

| Operation | Existing authority | Inputs → outputs and limits |
|---|---|---|
| nmr-deconvolution | Existing Wasserstein LP | Common uniform ppm grid, nonnegative spectra, explicit proton counts → known-component mole fractions/noise/distance ppm; dependent references rejected, no unknown identification |
| edi-transport | Existing parse_transport | Exact 5-column supplied table → T K and mobility cm²/V/s; no QE/EDI execution |
| edi-lifetimes | Existing parse_inv_tau/summarize | Exact 7-column supplied table → state indices, E eV, printed rate Ry, lifetime fs and summary |
| edi-matrix | Existing parse_edmat | Exact 9/13/15-column table → printed matrix records (M in Ry, squared M in Ry²); provider normalization retained, no guessed unit conversion |
| wannier-spreads | Existing parse_wout | Complete supported final text block → centers Å/spreads Å²; localization convergence requires review |
| epw-coupling-table | Existing parse_rows | Explicit selected 7-column table → band/mode indices, eV energies, meV frequency/coupling; no automatic CBM/gauge-invariance conclusion |
| uniform-kpoints | Existing uniform_kpoints | Positive integer mesh → full unshifted reciprocal-fraction grid; at most 10000 points, no convergence claim |
| xrd | pymatgen XRDCalculator | Ordered cell Å, wavelength Å, 2θ degrees → unbroadened peak positions, relative intensity, d Å, hkl; no Rietveld fit/profile broadening |
| structure-match | pymatgen StructureMatcher | Two ordered cells and explicit tolerances → equivalence; no database-wide novelty |
| supercell | pymatgen Structure.make_supercell | Cell and positive integer repeats → ordered cell; at most 2000 sites |
| substitute | pymatgen Structure.replace_species | Cell and full species map → ordered cell; no fractional occupancy/charged defect energetics |
| symmetry | pymatgen SpacegroupAnalyzer | Cell, length/angle tolerances → space group and conventional cell |
| slab | pymatgen SlabGenerator | Cell, hkl, thickness/vacuum Å, shift → unrelaxed slab and area Å²; polarity/termination not accepted automatically |
| hull | pymatgen PhaseDiagram | Uniquely identified phases with total eV and common basis → hull/formation eV per atom and decomposition; supplied phase completeness needs review |
| eos-fit | ASE EquationOfState | Distinct volumes Å³, total eV, common basis → Birch–Murnaghan V0/E0/B GPa; physical minimum must lie inside scan |
| ideal-gas | ASE IdealGasThermo | Geometry Å, complete positive vibrational energies eV, T K, P Pa, symmetry number, spin S → H/G eV and entropy eV/K; rigid-rotor harmonic approximation |
| rdf | ASE geometry.rdf.get_rdf | Periodic frames Å, atomic-number pair, cutoff Å, bins → normalized r Å and g(r); cell size enforced by ASE, no diffusion conclusion |
| probability-density | pymatgen ProbabilityDensityAnalysis | Fixed orthogonal-cell frames Å, species, spacing Å, explicit dt ps → raw normalized density Å⁻³; at most 50000 grid points, no smoothing/drift inference |
| intercalation-voltage | Existing skill calculate_voltage | Paired total energies eV and metal/ion counts → voltage V; restricted to monovalent Li/Na/K, stoichiometry/reference require review |
| grain-boundary-energy | Existing compute_gb_energy | GB eV, atom count, bulk eV/atom, actual cross-product area Å² → J/m² for exactly two equivalent interfaces |
| elastic-moduli | pymatgen ElasticTensor | Positive definite symmetric stiffness 6×6 GPa → VRH bulk/shear/Young GPa and Poisson ratio; no strain simulation |
| magnetic-moments | Existing magnetic parser | Species and collinear projected moments μB → site summary/sum/order heuristic; negative ferro sign fixed, no interstitial-inclusive moment or ground-state claim |
| dielectric-chi | Existing CHI parser | Unique frequency-dependent IPA/RPA block text → energy eV and printed scalar real/imaginary arrays; static blocks rejected, inverse-header convention retained |
| random-cubic-structure | Existing random generator | Integer composition, volume scale, seed → bounded distance-screened P1 candidate in cubic metric; no requested space-group symmetry or stability claim |
| molecular-descriptors | Existing RDKit utility | Valid SMILES and TPSA option → MW g/mol, TPSA Å², descriptors/heuristics; no clinical or experimental ADMET endpoint |
| molecular-fingerprints | Existing RDKit utility | SMILES, radius/bits/chirality/features → Morgan bits and Tanimoto matrix; source Butina bugs prevent exposing clustering |
| docking-box | Existing compute_box | Explicit selected coordinates Å, padding/minimum size Å → center/size Å; selection is the caller's responsibility |
| redocking-rmsd | Existing RDKit CalcRMS helper | Two bounded 3D molblocks with exact identity/stereochemistry → heavy-atom in-place RMSD Å; receptor frame preserved, no rigid alignment/PDBQT conversion |
| orca-energy | Existing parse_energy | Normally terminated text with decimal final energy → Hartree/eV; no executable launch, spin/orbital parsing or SCF convergence acceptance |
| error-metrics | Existing evaluate_metrics | Paired nonempty predictions/targets in stated eV/atom, eV/Å or eV/Å³ → MAE/RMSE; common label/reference basis explicit |
| spectrum-similarity | Existing spectrum functions | Nonflat spectra on a common increasing uniform ppm/cm⁻¹ axis → normalized l2/cosine/Wasserstein score; no identification/deconvolution |

The existing skill functions were moved into installed modules and re-exported by
source scripts/utilities, so formulas have one owner. Full tutorials may still need
unavailable providers. For example, the box-only operation does not install or
qualify MDAnalysis receptor selection. Unchanged pre-existing CLI issues are not
represented as successful end-to-end tutorial execution.

`identity(tool)` hashes the fixed provider source files and records the science,
jsonschema and relevant engine package versions. The recipe adds its own adapter
hash; each request is generated in its selected worker and must execute with the
same identity. The input artifact, specification, report and manifest remain linked
through the existing Job/attempt/receipt ledger. Scalar-looking outputs alone do
not justify creating canonical scientific results without registered targets.

Verification uses independent Bragg spacing, analytic EOS, ideal-gas and isotropic
elastic references, identified hull polymorphs, exact scalar/metric examples,
normalized occupancy, reproducible random candidates, and rigidly translated poses.
Invalid inputs, missing engines, implementation mismatch, malformed collection,
execution replay and separate site identities are covered. Worker tests skip RDKit
only in the numerical environment; a separate RDKit environment executes those
cases. See the integration receipt for exact counts and installed versions.

[environment-availability.json](environment-availability.json) records metadata
from 22 inspected environments, with private paths removed. A version entry is not
an import or execution certificate; null is absence in that environment, not a
claim about every possible install. Optional analysis/chemistry/transport extras are
isolated; package import and DB validation do not import heavy scientific engines.
Remote acquisition, licensed DFT/docking/phonon engines, generative checkpoints,
training pipelines and scientific acceptance decisions remain explicit prerequisites
or follow-up work in the complete matrix. They are not replaced with fake recipes.

Regenerate source inventory with `python tools/audit_skills.py` from this checkout
(PyYAML required). The checked-in review TSV is the semantic classification input;
the generator checks exact set equality against all 130 SKILL.md files and records
actual source hashes rather than inferring availability from skill names.
