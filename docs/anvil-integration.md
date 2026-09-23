# Optional external consumers

AtomisticSkills remains a standalone collection of workflows, skills, scripts
and MCP tools. Its normal research workflow does not require Anvil, a Campaign,
a database connection or a project distribution installation.

An external harness can read skill instructions and invoke existing scripts or
MCP tools. The consumer owns planning, authorization, durable execution state,
input bindings and result storage. Anvil-specific admission and publication rules
live in Anvil's `anvil_core.integrations.atomistic_skills`, not in this repository.

The Anvil adapter invokes `tools/run_analysis.py` in a selected scientific Python
environment. This optional JSON CLI calls the same helpers as the skill scripts;
it neither connects back to Anvil nor creates jobs. Its protocol contains explicit
scientific inputs and outputs, never database identities or campaign policy.
The integration audit is maintained by anvil-recipes under `docs/atomistic-skills`.

The previous `atomistic-analysis` distribution has been removed. Shared numerical
functions, timing validation, fit diagnostics and plotting remain source utilities
in `src/utils/analysis`. Existing skill script paths remain available. Numerical
dependencies still belong to the scientific environments documented by each skill.

Consumers upgrading from the old wheel must qualify a matching adapter and pinned
source revision. Historical outputs remain immutable; this source change does not
reinterpret earlier results or modify a running consumer environment.
