# MS/MS Agent Environment

This environment supports the original `coleygroup/ms-pred` tandem mass spectrum
models. It has no MCP server. The source environment name is `ms-gen`; the skill
catalog directory is `msms-agent`.

## Linux CPU installation

Create an isolated prefix from `core_env.yaml`, then activate that prefix. Run
these commands from this environment directory; set `MS_PRED_REPO_DIR` to the
provider checkout. The recorded provider revision is
`318bb7b743f67578115d6a4b5b2aae0ac62e3730`.

```sh
PYTHONNOUSERSITE=1 python -m pip install torch==2.4.0+cpu --index-url https://download.pytorch.org/whl/cpu
PYTHONNOUSERSITE=1 python -m pip install dgl --find-links https://data.dgl.ai/wheels/torch-2.4/repo.html -c provider-constraints.txt
PYTHONNOUSERSITE=1 python -m pip install torch-scatter torch-sparse --find-links https://data.pyg.org/whl/torch-2.4.0+cpu.html -c provider-constraints.txt
PYTHONNOUSERSITE=1 python -m pip install -e "$MS_PRED_REPO_DIR" -c provider-constraints.txt
PYTHONNOUSERSITE=1 python -m pip check
LD_PRELOAD="${CONDA_PREFIX:?}/lib/libstdc++.so.6" PYTHONNOUSERSITE=1 python -c 'import torch, dgl, ms_pred.common, ms_pred.nn_utils'
```

The provider needs RDKit 2025.03 and Ray 2.49.1 or newer; the former Ray 2.7 pin
was incompatible with this source revision. Keep the Torch constraint during
provider installation so the resolver cannot replace the declared 2.4 runtime
with a version incompatible with DGL/PyG extensions. Match CUDA wheels separately
when qualifying a GPU host. Import checks do not qualify trained models or
scientific spectrum predictions; follow the selected skill's checkpoint and
input instructions before advertising execution readiness.

On hosts with an older system C++ runtime, preload only the environment's
`libstdc++.so.6` as shown above when running model scripts. Adding the entire
Conda `lib` directory to `LD_LIBRARY_PATH` can instead select incompatible
Conda libtorch files ahead of the PyTorch wheel. Validate the imports together
in one process: separate import probes can miss that load-order conflict.
