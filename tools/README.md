# Source command line tools

`run_analysis.py` dispatches supplied-data operations to existing scientific
functions. It has no Anvil dependency; `--help` lists its interface.

`provision_environments.py` reads every environment's `core_env.yaml` (or
`env.yaml`) and installs into new explicit prefixes. Run with Python containing
PyYAML and an existing Conda installation:

```sh
python tools/provision_environments.py plan --directory "$INSTALL_ROOT"
python tools/provision_environments.py install --directory "$INSTALL_ROOT" \
  --conda "$CONDA_EXECUTABLE" --workers 2
```

Use `--environments base-agent phasefield-agent` to select a subset. Existing
environment directories are refused, including previous failed attempts. Choose
a fresh directory for retries; previous logs remain evidence. This tool never
removes a named environment or runs the legacy `install.sh` scripts. A failed
environment makes the command exit nonzero; inspect each `receipt.json` and log.

Dependency installation is not engine or scientific qualification. Follow each
environment README for extra provider checkouts, extension wheels, checkpoints,
licensed binaries and device requirements. Record those steps and run actual
tool tests before advertising a skill as executable. `pip check` alone cannot
validate CUDA drivers, model compatibility, external services or scientific
outputs. MCP servers currently import the MCP 1.x SDK; the declarations bound
that major version to prevent the incompatible 2.x module layout.
Generated MCP configuration sets `PYTHONNOUSERSITE=1` so unrelated packages in the
user's home directory cannot override the chosen Conda environment. For manual
script calls set the same variable. External provider paths use `ADIT_REPO_DIR`,
`DIFFCSP_REPO_DIR` and optionally `DIFFCSP_CHECKPOINT_DIR`; atomate2 requires an
explicit `ATOMATE2_VASP_CMD` or `VASP_CMD` instead of guessing a developer's path.

The provisioning tests preserve existing prefixes and exercise resolver failures.
`tests/base/test_mcp_json_arguments.py` starts a fresh real server and checks JSON
string and structured matrix arguments; run it in `base-agent`. Restart existing
MCP server processes after source changes to load the corrected signatures.
