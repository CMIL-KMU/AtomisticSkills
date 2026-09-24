# Peregrine MACE fine-tuning

Use a qualified `peregrine-mace` environment containing Peregrine, MACE, PyTorch,
Lightning, ASE and pymatgen. The script imports these installed packages; it does
not install packages or download checkpoints. Supply a trusted local MACE model.

```bash
# Env: peregrine-mace
python .agents/skills/ml-mace-finetune/scripts/peregrine_finetune.py \
  --dataset dataset.json --foundation foundation.model --output-dir results \
  --seed 42 --epochs 100 --lr 0.0001 --batch-size 2 --device cuda
```

The output directory must be empty. `--trainable readout` (default) freezes the
backbone, preserving foundation readouts as initialization; `--trainable all`
updates all parameters. `--head NAME` selects a foundation head where supported.
`--provider-sha256 HASH` additionally rejects changes to installed Peregrine Python
source files. The report always records those hashes.
Allocate a GPU through your site's scheduler before selecting CUDA.

Input JSON has `units` exactly equal to
`{"energy":"eV","forces":"eV/angstrom","stress":"eV/angstrom3","stress_sign":"tensile"}`
and a `structures` list. Each entry contains a unique `id`, an explicit `split`
(`train`, `validation`, or `test`), a pymatgen `structure` dictionary, scalar total
`energy`, N×3 `forces`, and symmetric 3×3 `stress`. All labels must be finite.
All splits must be nonempty; exact duplicate geometries are rejected. The caller
must additionally group related configurations to avoid scientific data leakage.
Raw VASP kbar stress needs sign reversal and division by 1602.1766208 upstream;
already canonical stress must not be converted again.

The held-out test frames never enter the trainer. The train/validation rows are
reordered to preserve the requested membership with Peregrine's native seeded
split. Foundation energy references are retained (`native_passthrough`, no
scale/shift fitting). Training uses float32, Adam, linear warmup/cosine decay,
and early stopping with patience 20. Losses are MSE of per-atom energy (weight 1),
force components (weight 10), and stress converted to GPa (weight 1).

`report.json` reports MAE/RMSE for each split: meV/atom energy, eV/Å forces and
GPa stress (six independent components). These are absolute energy errors without
an offset fitted on test labels. Improvement is not guaranteed. `predictions.json`
preserves targets and before/after predictions. The selected model is exported as
`finetuned.model` and `weights.npz`; the native export is reloaded to check E/F/stress
agreement. Input/output hashes, installed provider source hashes, package versions,
training configuration and CSV history are retained. Retain the full output folder
for provenance. A small correlated dataset only verifies this workflow, not
transferability or response outside its sampled structures and strains.
