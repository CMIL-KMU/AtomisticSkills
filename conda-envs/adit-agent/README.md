# ADiT Agent Environment

Environment for the All-atom Diffusion Transformer (ADiT) — a unified latent diffusion model for generating both periodic crystals and non-periodic molecules.

## Prerequisites

- Conda/miniforge3 installed
- Python 3.10
- CUDA 13.0 compatible GPU

## Installation

1. Create the conda environment:
```bash
conda env create -f core_env.yaml
conda activate adit-agent
```

2. Install PyTorch 2.9.1+cu130:
```bash
pip install torch==2.9.1 torchvision --index-url https://download.pytorch.org/whl/cu130
```

3. Install PyG and extensions:
```bash
pip install torch-geometric
```

> **GB10 / ARM / CUDA 13.0 users**: torch-scatter and torch-cluster must be compiled from source. Follow the **"Installation on ARM"** section in [`conda-envs/mattergen-agent/README.md`](../mattergen-agent/README.md) for `CUDA_HOME`, `TORCH_CUDA_ARCH_LIST`, and build commands.

4. Install remaining dependencies:
```bash
pip install lightning==2.4.0 hydra-core hydra-colorlog
pip install e3nn==0.5.1 einops rootutils rich omegaconf torchdiffeq huggingface_hub
pip install pymatgen ase rdkit pyxtal tqdm scipy pandas matplotlib torchmetrics
pip install timm lmdb wandb pathos p-tqdm download
pip install smact matminer importlib_resources posebusters
pip install svgwrite CairoSVG reportlab svglib pythreejs ipywidgets
pip install mofchecker pymatgen-analysis-defects orjson lovely-tensors submitit
conda install -c conda-forge openbabel
pip install 'mcp>=2.2,<3'
```

5. Clone the AADT repository:
```bash
git clone https://github.com/facebookresearch/all-atom-diffusion-transformer "$ADIT_REPO_DIR"
```

## MCP Tool

| Tool | Description |
|---|---|
| `generate_structures` | Generate crystals (MP20) or molecules (QM9) unconditionally |

**Parameters**: `generation_type` (`crystals` / `molecules`), `num_structures`, `cfg_scale`, `batch_size`, `device`, `output_dir`

## Checkpoints

Pre-trained checkpoints are auto-downloaded from HuggingFace on first use:
- Repository: `chaitjo/all-atom-diffusion-transformer`
- `ckpts/ldm.ckpt` — Latent diffusion model
- `ckpts/vae.ckpt` — Variational autoencoder

## Key Dependencies

| Package | Version |
|---|---|
| Python | 3.10 |
| torch | 2.9.1+cu130 |
| torch-geometric | 2.7.0 |
| lightning | 2.4.0 |
| e3nn | 0.5.1 |
| pymatgen, ase | latest |

## Provider validation

Install the evaluation dependencies above even for inference: upstream imports
its evaluators while loading the model. This avoids modifying scientific provider
code or silently disabling evaluators. Set `ADIT_REPO_DIR` to the cloned repository
and `HF_HOME` to an explicit model cache. Test the original wrapper from the
AtomisticSkills repository root:

```bash
PYTHONNOUSERSITE=1 python -c 'from src.utils.generative_models.adit.adit_wrapper import ADiTWrapper; ADiTWrapper(device="cpu")'
```

CPU model initialization verifies imports and checkpoint compatibility; it does
not qualify generation speed or CUDA execution. Use the intended GPU for the
actual generation qualification.

## Notes

- Set `ADIT_REPO_DIR` to an explicit checkout path, or use the sibling `adit` discovery convention.
- First run downloads ~3 GB of model checkpoints from HuggingFace
- ADiT supports **dataset-type** (crystals vs molecules) and **spacegroup** conditioning only — no composition or property conditioning
