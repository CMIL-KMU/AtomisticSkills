# Independent simulation source tools

`tools/run_simulation.py ACTION --input request.json --output response.json` accepts JSON
and produces one JSON response. Omitting the file flags uses stdin/stdout. Diagnostics
from Python/native providers go to stderr. Output files and calculation directories must
be new. No project distribution or workflow framework is required.

## Environments

Peregrine uses the engine's own locked environment installer. The qualified integration
revision is `66c66326046397ea7560f680042eb6078f43369a` from `CMIL-KMU/peregrine`.
In that checkout, use `./env.sh install core-cpu` for CPU or the engine's documented UMA
profile with its explicitly licensed checkpoint. Install this source tool's `jsonschema`
validation dependency in the same isolated prefix. Do not install the legacy Peregrine
agent. GPU execution requires engine hardware qualification in the selected environment.

## Peregrine

`peregrine-run` accepts `{ "request": {...}, "registry": "/absolute/model-registry",
"output": "/absolute/new-output" }`. The physical request schema is `peregrine/request.json`.
It has no job IDs, DB references, or scheduler settings. `device` defaults to `cpu`.

Results include `simulation.json`, a staged `deployment` registry, and applicable
checkpoint, trajectory, frequency, path, selection or bias artifacts. The response
retains energy, force and stress conventions, model identity, outcome, constraints,
ordered system IDs, source bytes and dependency versions. These are scientific
observations; only a consuming workflow can accept them into its DB contract.

New foundation descriptors use `peregrine.uma-foundation/v1`. Legacy hash-bound
descriptors remain readable for saved-model identity compatibility.

## Qualification

`tests/simulations` checks XML observations and restricted foundation loading.
For the real standalone CLI test, set `ATOMISTIC_PEREGRINE_FIXTURE` to a directory
containing an engine-created `registry/` and `fixture.json` with `model`, ordered
`structures`, and reference `expected` energy/force/stress rows. The test blocks all
consumer workflow imports and compares real provider calculations to those reference rows.
Keep test output under `.agents/test` using pytest's `--basetemp` option.
