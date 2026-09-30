# Scientific source utilities

These functions are shared by existing skill scripts and the optional
`tools/run_analysis.py` JSON CLI. They require no project wheel, database,
scheduler or agent framework. Run scripts in the environment specified by the
skill. Direct Python consumers can use `src.utils.analysis` from this checkout.

The original skill entry points remain:

```sh
python .agents/skills/mat-diffusion-analysis/scripts/analyze_diffusion.py --help
python .agents/skills/mat-diffusion-analysis/scripts/calculate_activation_energy.py --help
echo '{}' | python tools/run_analysis.py describe
echo '{"tool":"uniform-kpoints","inputs":{"mesh":[2,2,2]}}' | python tools/run_analysis.py run
```

The CLI dispatches only named actions. `describe` exposes scientific input schemas;
`identity` records source hashes and numerical dependency versions. It does not
admit operations for a consumer or publish results into any database.

Diffusion retains explicit species, charge, temperature and timing, fixed-cell
validation and single-origin MSD defaults. Contradictory timing and unsupported
max-smoothed frame intervals are rejected. Arrhenius retains exclusions, regression
diagnostics and opt-in extrapolation. Numerical success does not imply scientific
acceptance. Existing outputs are not rewritten.

Transport plotting requires pymatgen-analysis-diffusion 2025.11.15, ASE,
Matplotlib, PyYAML and LovelyPlots 1.0.2 in the chosen environment. Structure PNGs
also require pymatviz 0.18.0, Plotly 7.0.0, Kaleido 1.4.0 and an explicitly
configured browser (`BROWSER_PATH`). `structure-png` accepts saved periodic
coordinates on stdin and writes PNG bytes; other actions return JSON.

`transport-contract` exposes defaults and parameter semantics; `transport-review`
resolves a supplied configuration and reports explicit inputs, inherited defaults
and deviations. Its digest binds effective settings to the scientific source.
Consumers decide whether a deviation needs review; this utility never approves it.

With default `smoothed=false`, plots show the full available original MD time:
a 0–30 ps trajectory with `equilibration_ps=5` displays 0–30 ps, referenced to
the first frame. The first 5 ps are shaded and excluded from the fit; they are
not deleted from the MSD curve. Fit bounds remain lag times from that reference
and are intersected with the post-equilibration interval. Multiple-origin
`smoothed="max"` still uses post-equilibration frames and lag-time plotting.
`msd.csv` keeps `lag_ps,msd_angstrom2,plot_time_ps` and appends `fit_included`.
The linear fit has a free intercept, recorded separately in the result.
`input_configs.yaml` records requested/effective settings, source contract, input
timing and the resolved analysis, plot and fit intervals. No old output is changed.

This full-time single-origin convention uses the first trajectory frame as its
reference, unlike the previous post-equilibration reference. It can change fitted
values; preserve the source identity and do not relabel earlier analyses.
