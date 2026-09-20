# Verification receipt, 2026-09-21

Branch `feat/skill-integration-sites`, science package version 0.2.0. Source pin is
recorded by the consuming recipe release after this source commit. No production
runtime or scientific Job was changed. Full private logs, environment metadata and
SSH observations remain in `.agents/test/integration-20260921`.

- Complete source census: 130 skills, 223 Python script contracts, 76 explicit MCP
  references; three unresolved references are marked in the inventory.
- Registry: 31 operations mapped to 29 partially integrated skills. Four other
  skills retain existing integrations. The remaining classifications are 61
  prerequisite, 14 acquisition, 9 guidance, 7 literature, 3 review-assisted,
  2 science-blocked and 1 reference-blocked. Partial never means full tutorial.
- Numerical/transport/Arrhenius/presentation regression: 43 passed, 2 RDKit skips
  in the installed Python 3.12 numerical environment. These counts include the
  provider subset (20 passed, 2 skips); do not add the overlapping counts.
- Chemistry: 3 passed in Python 3.12 with separately installed RDKit 2025.9.6.
  An earlier Python 3.9 / RDKit 2024.3.3 exploratory pass is not the supported
  package-runtime qualification. The box test overlaps the numerical suite.
- Existing recipe transport tests: 4 passed. Native installed recipe plus
  synthetic scheduler/workflow checks: 66 passed. Additional venv analysis-profile
  execution: 1 passed. Existing CPU/GPU profiles and registry: 3 passed in a real
  Conda environment; attempting these Conda-only profiles in a venv correctly
  failed before the compatible analysis profile was added.
- Source lint on new modules/tests and whitespace checks passed. The wider legacy
  source tree has existing compact-style lint and unrelated CLI issues; no claim
  of a clean whole-repository Ruff run is made.

Tests cover independent Bragg/EOS/hull/ideal-gas/tensor/scalar references, occupancy
normalization, known NMR mixture, exact EDI/EPW fields, pose translation, invalid
input, absent engine, source hash mismatch, output replay and missing worker.
Pymatgen core's split distribution is recorded independently; older releases report
that core is bundled in the identified pymatgen distribution. Heavy engines remain
lazy. The consumer's DB acceptance and full workspace gate receipts live in Anvil's
`spikes/skill-integration/RECORD.md`; they are separate from numerical qualification.
