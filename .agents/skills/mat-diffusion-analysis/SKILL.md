---
name: mat-diffusion-analysis
description: Analyze ionic self diffusion, Nernst–Einstein conductivity and Arrhenius activation energies from explicit MD trajectories or observations.
category: [materials]
---

# Diffusion Analysis

## Goal
Use the shared `src.utils.analysis` tools to obtain reproducible provisional
transport estimates and identify evidence missing for scientific acceptance.

## Instructions

1. Identify the source MD and its saved-frame spacing, equilibrated interval,
   mobile species, ionic charge and temperature. Reuse the approved research
   settings and ask only for missing scientific choices. Select existing MD tools
   from the relevant simulation skill when a new trajectory is needed.
2. Use the existing skill scripts or source utility API. Default `smoothed=false` uses one time
   origin; `max` uses multiple origins. Select a diffusive fit interval after
   inspecting MSD. Short trajectories and apparent straight lines do not establish
   convergence. Do not treat a live monitor's early stop as validated transport.

```bash
# Env: base-agent
python .agents/skills/mat-diffusion-analysis/scripts/analyze_diffusion.py trajectory.traj --species Li --charge 1 --temperature 600 \
  --frame-interval-fs 10 --ignore_ps 5 --smoothed false --output_dir analysis-600K
```

For physical_step-tagged Peregrine extxyz, use `--timestep-ps 0.001` instead of
`--frame-interval-fs`. Fit lag limits: `--fit-start-ps` and `--fit-end-ps`. CLI output
must be a new directory. Effective/default settings are in input_configs.yaml;
JSON, MSD CSV and PNG/SVG preserve the evidence.

3. Select compatible temperature observations explicitly. Each JSON row contains
   source_uuid, temperature_K, diffusivity_cm2_s and optional conductivity_NE_mS_cm.
   Preserve species, ionic_charge_e and method if available; mixed values are rejected.
   Retain source identities; never infer comparability from folder names.

```bash
# Env: base-agent
python .agents/skills/mat-diffusion-analysis/scripts/calculate_activation_energy.py --observations observations.json --output_dir arrhenius-analysis
```

4. Inspect exclusions, eligible temperatures and residuals. Nonpositive, missing or
   nonfinite coefficients are excluded with reasons, not deleted from their source.
   Repeated temperatures need explicit replicate handling. Three eligible temperatures
   are required for regression scatter. Nernst–Einstein fitting uses sigma*T.
5. Ask whether prediction at another temperature is wanted and which temperature.
   Only then pass `--target-temperature-K 298.15` (example). There is no fixed Li
   species or automatic 300 K extrapolation. Predictions assume the same mechanism.
   Retain source hashes and method identity with every reported result.

## Examples

The upstream [LGPS study](examples/LGPS/README.md) is historical evidence. Its
numbers used the upstream environment/protocol; do not overwrite or present them
as validation of a new method. See the [API contract](../../../src/utils/analysis/README.md)
for explicit input fields and migration from the previous CLI.

## Constraints

Fixed-cell periodic solids with mobile ions and a nonempty framework are supported.
Atom order must be preserved and frame timing explicit/uniform. NPT is not currently
qualified. Resolve ion charge explicitly; Nernst–Einstein ignores cross correlations.
The tool rejects nonpositive raw fitted slopes before pymatgen's positive floor.
Fit standard errors are regression diagnostics, not independent-sample uncertainty.
All scripts call the common source utilities; do not implement another estimator in a skill.

## References

- pymatgen-analysis-diffusion, DiffusionAnalyzer and get_diffusivity_from_msd.
- [Integration responsibilities](../../../docs/anvil-integration.md).

---

**Author:** Bowen Deng (upstream); laboratory integration requested by Hoje Chun.
**Contact:** [GitHub @CMIL-KMU](https://github.com/CMIL-KMU)
