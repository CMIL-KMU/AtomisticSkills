# Independent simulation source tools

`tools/run_simulation.py ACTION --input request.json --output response.json` accepts JSON
and produces one JSON response. Omitting the file flags uses stdin/stdout. Diagnostics
from Python/native providers go to stderr. Output files and calculation directories must
be new. There is no AtomisticSkills wheel or Anvil/mkite import requirement.

## Environments

VASP preparation/parsing uses `conda-envs/base-agent/core_env.yaml` (ASE, pymatgen,
Jinja2 and PyYAML). The caller supplies licensed PAW files, a VASP executable and any
scheduler. This source tool does not submit jobs or decide computational retries.

Peregrine uses the engine's own locked environment installer. The qualified integration
revision is `66c66326046397ea7560f680042eb6078f43369a` from `CMIL-KMU/peregrine`.
In that checkout, use `./env.sh install core-cpu` for CPU or the engine's documented UMA
profile with its explicitly licensed checkpoint. Install this source tool's `jsonschema`
validation dependency in the same isolated prefix. Do not install the legacy Peregrine
agent. GPU execution requires engine hardware qualification in the selected environment.

## VASP

- `vasp-registry`: `{ "protocol_directory": "/absolute/catalog" }` returns registry and aliases.
- `vasp-prepare`: provide `protocol_directory`, `name`, a pymatgen `structure` dictionary,
  `overrides`, explicit ambiguity `selections`, and optional `dynamics`/`ldau_preset`.
  Additional dynamics geometry bindings contain only `{ "structure": ... }`.
- `vasp-inputs`: provide `folder`, `structure`, prepared `stage`, ordered `potcar_paths`,
  optional machine `parallel` INCAR tags and `previous` directory for CHGCAR continuation.
- `vasp-parse`: provide `directory` for a complete single-stage output; provide `directory`
  and `stage` for dynamics-aware parsing (see `vasp/dynamics.py`).
- `identity`: `{ "provider": "vasp" }` records implementation hashes/dependencies.

Protocol catalogs remain caller-supplied until redistribution rights are established.
No licensed PAW datasets or VASP binaries are included. Existing profile aliases are
accepted for saved-input compatibility. The standalone DFT parser and this interface
share `ionic_steps.e_wo_entrp` energies and tensile-positive stress in eV/angstrom^3.
OUTCAR supplements run statistics; it does not overwrite XML force/stress observations.

## Peregrine

`peregrine-run` accepts `{ "request": {...}, "registry": "/absolute/model-registry",
"output": "/absolute/new-output" }`. The physical request schema is `peregrine/request.json`.
It has no job IDs, DB references, or scheduler settings. `device` defaults to `cpu`.

Results include `simulation.json`, a staged `deployment` registry, and applicable
checkpoint, trajectory, frequency, path, selection or bias artifacts. The response
retains energy, force and stress conventions, model identity, outcome, constraints,
ordered system IDs, source bytes and dependency versions. These are scientific
observations; only a consuming workflow can accept them into its DB contract.

Legacy hash-bound `anvil.uma-foundation/v1` model descriptors remain readable without
an Anvil dependency. New descriptors may use `peregrine.uma-foundation/v1`.

## Qualification

`tests/simulations` checks XML observations and restricted foundation loading.
For the real standalone CLI test, set `ATOMISTIC_PEREGRINE_FIXTURE` to a directory
containing an engine-created `registry/` and `fixture.json` with `model`, ordered
`structures`, and reference `expected` energy/force/stress rows. The test blocks all
Anvil/mkite imports and compares real provider calculations to those reference rows.
Keep test output under `.agents/test` using pytest's `--basetemp` option.
