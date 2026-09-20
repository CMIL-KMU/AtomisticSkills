---
trigger: always_on
description: Scientific protocols and Anvil campaign boundaries
---

# Research standards for the laboratory fork

AtomisticSkills owns scientific methods, skill instructions and reusable tools.
Anvil owns Project, Campaign, Iteration, requester, approval/advancement policy,
Job/result provenance and research state. anvil-recipes binds versioned tools to
Anvil inputs/outputs; anvil-machine owns cluster submission and monitoring.

For questions and reviews, answer directly. For a computational research task:

1. Reuse the user's existing Campaign and its approved settings when provided.
   Consult protocol documents for scientific guidance; do not create a parallel
   campaign in local Markdown/CSV or assume authority to submit jobs.
2. Propose missing scientific settings and record the approved plan through the
   available Anvil integration. If the integration is unavailable, return a plan
   proposal; do not invent an API or bypass it with direct scheduler submission.
3. Run approved steps through recipes/machine. Use the allocated Job workspace
   and retain explicit source/result identities. Local tools can also be called
   directly for standalone analyses when the user requests that scope.
4. Judge diffusion regime, fit interval, phase stability and model validity.
   Report incomplete evidence as provisional; numerical success is not validation.

The documents in `.agents/workflows/` are scientific protocol references, not an
independent execution state machine. Campaign state and human review decisions
live in Anvil. Markdown reports and CSV exports are views, not a second authority.
Do not require a separate research_plan.md approval after the user has already
approved the same work. Literature lookup and model selection remain scientific
activities, performed when relevant rather than for every computation.
