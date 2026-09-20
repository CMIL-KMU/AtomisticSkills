# Laboratory integration boundaries

Anvil constructs Campaigns and owns Project/Iteration state, request attribution,
approvals, results, datasets and lineage. anvil-recipes translates versioned tool
inputs/outputs. anvil-machine submits and monitors cluster work. Peregrine supplies
MLIP models and torch-sim simulation. AtomisticSkills owns reusable scientific
protocols, interpretation guidance and analysis tools.

The upstream `.agents/workflows` documents remain scientific references. They
must not create a competing CSV/Markdown campaign ledger. The laboratory rules
replace the mandatory standalone research-plan/approval ceremony with the current
Anvil Campaign policy. Local analysis remains possible outside Anvil when requested.
Existing atomate2/jobflow-remote tools are retained for upstream compatibility;
they are not the laboratory Campaign submission route. Do not introduce another
scheduler monitor or duplicate Peregrine model wrappers for these Campaigns.

Transport and Arrhenius now share the installable `atomistic_analysis` package.
The skill CLI and Anvil recipes are adapters to this package. Managed analysis
uses explicit source Job/attempt/hash inputs, never directory discovery. Normal
CLI use accepts a trajectory or explicit observation JSON, without a DB.

New transport analyses default to single-origin (`smoothed=False`). Historical
max-smoothed results remain immutable. Changing method/fit range requires a new
analysis revision/Job, and mixed methods cannot be silently combined. Excluded
coefficients and finite-trajectory uncertainty remain visible. No scientific
acceptance is implied by API success. Deployment must install the same versioned
analysis wheel in worker and reader environments; Arrhenius readers need only the
lightweight base package. Preserve commit/wheel hashes in the release receipt.

Fork development: keep `upstream` pointing to learningmatter-mit/AtomisticSkills,
`origin` to the laboratory fork, and feature branches in separate worktrees.
Preserve upstream attribution and the MIT license. Do not commit research inputs,
trajectories, credentials or generated environment directories.
