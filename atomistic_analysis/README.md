# Scientific analysis API

This installable package owns deterministic transport and Arrhenius calculation.
It has no Anvil, MCP, MLIP or scheduler dependency. Arrhenius uses only the standard
library; trajectory analysis uses the optional `transport` extra.

In the existing base-agent Conda environment, from this repository root:

```sh
python -m pip install '.[transport]'
atomistic-diffusion --help
atomistic-arrhenius --help
```

Python: `atomistic_analysis.transport.analyze(path, config, temperature_K=600,
frame_interval_fs=10)`; configuration requires species, ionic charge and discarded
equilibration time. Physical-step trajectories use `timestep_ps` instead. Do not
supply both timing conventions. `smoothed=False` is the default; `"max"` remains
available. `fit_start_ps`/`fit_end_ps` select lag times after equilibration removal.
The single-origin reference is the first retained frame. The default fit uses
all lag points from the selected estimator; inspect and restrict to the diffusive
regime when necessary. MSD and regression standard error do not establish convergence.

`atomistic_analysis.arrhenius.fit_transport(rows, target_temperature_K=None)`
accepts explicitly selected source identities, temperature_K, diffusivity_cm2_s
and optional conductivity_NE_mS_cm. It excludes missing/nonpositive/nonfinite
coefficients with reasons, rejects duplicate trajectory IDs and repeated eligible
temperatures, and needs three eligible temperatures. D and sigma*T are fitted in
log10 space with unweighted least squares. Optional target temperature returns
log10 predicted values and their units without overflow or implicit 300 K selection.
No extrapolation occurs by default. These are Nernst–Einstein estimates, not
collective charge conductivity; ion-ion cross correlations are not modeled.

One source of numerical truth: skill CLIs and anvil-recipes call this API.
Anvil result schemas are applied only by the recipes adapter. API results identify
the provider version, effective configuration and numerical implementation hash.
CLI effective inputs are saved separately in input_configs.yaml. Earlier CLI
folder-scanning and implicit timing/charge/300 K behavior is intentionally replaced
with explicit inputs; old results remain historical data, not reinterpreted.

The initial supported trajectory contract is fixed-cell periodic ionic solids
with a nonempty framework and preserved atom order. NPT and pure single-species
liquids need separate validation; they are rejected here. Periodic unwrapping
requires sufficiently frequent frames; a full cell traversal cannot be detected
from wrapped frames alone.

## Presentation

Scientific exports use `plotting.plot_style()` from LovelyPlots 1.0.2 (`ipynb`).
`plotting.browser_style()` exposes the same palette and style metadata for browser
charts. The installed upstream stylesheet is read through Matplotlib's public file
API because LovelyPlots' registration hook uses a removed private Matplotlib API.
Historical saved plots are immutable; new renders use the shared style.

`structure_viz.structure_3d_custom` is the original AtomisticSkills Plotly renderer,
now importable from the wheel; the old utility path re-exports it. Install the
`visualization` extra and an explicitly configured Chrome/Chromium executable
(`BROWSER_PATH`) for PNG generation. `python -m atomistic_analysis.structure_cli`
accepts species, fractional coords and lattice JSON on stdin and writes PNG to
stdout. It displays saved periodic coordinates, no inferred bonds, in Cartesian
x/y/z and perspective views; no structure standardization occurs. Preview limit
is 500 sites. Nonperiodic molecule rendering is outside this endpoint's scope.
