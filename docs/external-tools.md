# External scientific tools

Agents can discover `.agents/skills/*/SKILL.md`, load relevant rules/workflows,
and invoke the existing scripts or MCP tools in each skill's environment. The
repository remains usable directly for interactive research. No project wheel,
database or external workflow service is required.

`tools/run_analysis.py` provides a JSON interface to shared numerical utilities
under `src/utils/analysis`. `tools/run_simulation.py` exposes the independent
Peregrine interface documented in `src/simulations/README.md`. VASP's existing
local scripts use `src/utils/dft/vasp_writer.py` and `vasp_parser.py`; reusable
strict stage extraction is available in `vasp_results.py`.

These interfaces contain scientific inputs, physical controls, observations,
units and provenance. A consuming harness owns durable jobs, scheduling,
database identity, acceptance policy and the next research decision. Consumers
should pin source revisions and qualify their input/output adapters. Discovery
of a skill does not by itself install its programs or establish execution support.

Common numerical corrections belong in these reusable helpers. Consumer-specific
adapters and unpublished protocol catalogs belong in their owning repositories.
Keep raw outputs and source identities so results can be independently inspected;
do not reinterpret saved results after a provider upgrade.
