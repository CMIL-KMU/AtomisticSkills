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
