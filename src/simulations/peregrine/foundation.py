"""Portable, hash-bound UMA foundation input; no caller-selected Python classes."""
import hashlib
import json
import math
import re
from pathlib import Path
import shutil


def load_foundation(root, identity, device):
    import torch
    from peregrine.artifacts.store import canonical, file_digest
    from peregrine.train.factory import build_model_from_spec
    from peregrine.train.unified_config import ModelSpec, BackboneSpec, ScaleShiftSpec
    from ase.data import chemical_symbols
    directory = Path(root) / identity
    spec = json.loads((directory / 'foundation.json').read_text())
    if set(spec) != {'schema', 'sha256', 'species', 'cutoff', 'task'}:
        raise ValueError('Unexpected foundation specification fields')
    if spec['schema'] not in {'peregrine.uma-foundation/v1', 'anvil.uma-foundation/v1'} or hashlib.sha256(canonical(spec)).hexdigest() != identity:
        raise ValueError('Foundation specification differs from model identity')
    if spec['task'] not in ('omat', 'oc20', 'oc22', 'omol') or not isinstance(spec['species'], list) or not spec['species']:
        raise ValueError('Explicit supported UMA task and species required')
    if len(set(spec['species'])) != len(spec['species']) or any(s not in chemical_symbols[1:] for s in spec['species']) or type(spec['cutoff']) not in (int, float) or not math.isfinite(spec['cutoff']) or spec['cutoff'] <= 0 or not isinstance(spec['sha256'], str) or not re.fullmatch('[0-9a-f]{64}', spec['sha256']):
        raise ValueError('Foundation species, cutoff and checkpoint digest must be explicit and valid')
    checkpoint = directory / 'checkpoint.pt'
    if checkpoint.is_symlink() or file_digest(checkpoint) != spec['sha256']:
        raise ValueError('Foundation checkpoint hash differs from portable specification')
    # The provider/vendor subsequently opens this same private checkpoint. Prove
    # it is a restricted checkpoint first; never allow an arbitrary pickle global.
    from collections import defaultdict
    from typing import Any
    from omegaconf import DictConfig, ListConfig, OmegaConf
    from omegaconf.nodes import AnyNode
    from omegaconf.base import ContainerMetadata, Metadata
    from fairchem.core.units.mlip_unit.api.inference import MLIPInferenceCheckpoint
    with torch.serialization.safe_globals([MLIPInferenceCheckpoint, DictConfig, ListConfig, AnyNode, ContainerMetadata, Metadata, dict, list, int, Any, defaultdict]):
        saved = torch.load(checkpoint, map_location='cpu', weights_only=True, mmap=True)
    if not isinstance(saved, MLIPInferenceCheckpoint):
        raise ValueError('Expected an audited MLIP inference checkpoint')
    targets = {'fairchem.core.models.base.HydraModel', 'fairchem.core.models.uma.escn_moe.eSCNMDMoeBackbone',
               'fairchem.core.models.uma.escn_md.MLP_EFS_Head', 'fairchem.core.models.uma.escn_moe.DatasetSpecificMoEWrapper',
               'fairchem.core.units.mlip_unit.mlip_unit.Task', 'fairchem.core.modules.loss.DDPMTLoss',
               'fairchem.core.modules.loss.PerAtomMAELoss', 'fairchem.core.modules.loss.L2NormLoss',
               'fairchem.core.modules.loss.MAELoss',
               'fairchem.core.modules.normalization.normalizer.Normalizer',
               'fairchem.core.modules.normalization.element_references.ElementReferences', 'torch.DoubleTensor'}
    def check(value):
        if OmegaConf.is_config(value):
            value = OmegaConf.to_container(value, resolve=False)
        if isinstance(value, dict):
            for key, child in value.items():
                if key in ('_target_', 'model', 'head_cls', 'module') and child not in targets:
                    raise ValueError('Unreviewed foundation reconstruction target')
                check(child)
        elif isinstance(value, (tuple, list)):
            for child in value:
                check(child)
        elif isinstance(value, str) and '${' in value:
            raise ValueError('Foundation configuration interpolation is not allowed')
    check(saved.model_config)
    check(saved.tasks_config)
    model_spec = ModelSpec(strategy='native_passthrough', scale_shift=ScaleShiftSpec(source='none'),
                          backbone=BackboneSpec(type='uma', source='pretrained', name_or_path=str(checkpoint), species=spec['species'], cutoff=spec['cutoff'], kwargs=dict(task_name=spec['task'], device=device)))
    potential = build_model_from_spec(model_spec, dtype=torch.float32).eval()
    potential.output_properties = ('energy', 'forces')
    return potential, dict(target=dict(selectors=dict(head=[0], fidelity_idx=[0])), model_identity=dict(model=identity, checkpoint_sha256=spec['sha256'], task=spec['task']))


def copy_foundation(source, identity, destination):
    shutil.copytree(Path(source) / identity, Path(destination) / identity)
