---
name: mat-mlip-simulation
description: Run reproducible periodic ML potential simulations with Peregrine, preserving ordered structures, explicit physical settings, numerical evidence and restart checkpoints.
category: materials
---

# Periodic ML potential simulations

Use an installed Peregrine pot/sim environment and a content-addressed model registry.
This skill runs independently through the source CLI; no database account or workflow service is required.

1. Read [simulation inputs and environments](../../../src/simulations/README.md).
   Select a qualified model and environment. A random-initialized test model is not a research potential.
2. Write a JSON request matching [the physical request schema](../../../src/simulations/peregrine/request.json).
   Set `operation`, ordered `structures`, `options`, `seed`, and the immutable `model` identity.
   Supported operations are static evaluation, fixed-cell relaxation, MD, NEB, frequencies,
   and sampling. Preserve species, atom order, constraints, selectors and units.
3. Wrap the request with explicit absolute `registry` and new `output` directory paths.
   Run from the AtomisticSkills checkout:

   ```bash
   # Env: qualified Peregrine pot/sim runtime; see src/simulations/README.md
   python tools/run_simulation.py peregrine-run --input invocation.json --output response.json
   ```

4. Read `simulation.json`, termination/convergence, units, source/dependency identity and
   operation-specific artifacts before interpreting results. ML predictions are not DFT labels.
   Preserve the complete output directory. A failed command or partial directory is not a completed result.
5. For MD/NEB continuation, point `registry` at the prior output's `deployment` directory,
   use the saved final geometries and `restart.manifest_sha256` from its checkpoint receipt.
   Preserve source/dependency identity, model, seed, system order and physical configuration.
   Choose a new output directory; do not overwrite an existing attempt.

An external harness may submit these same inputs to a workflow adapter and read its validated
results before deciding the next step. The harness owns research decisions; this skill owns
physical input/output behavior and reuses the installed Peregrine engine.

---

**Author:** Hoje Chun (laboratory integration).
**Contact:** [GitHub @CMIL-KMU](https://github.com/CMIL-KMU)
