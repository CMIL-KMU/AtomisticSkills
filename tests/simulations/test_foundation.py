"""Restricted UMA deserialization counterexamples in an explicitly installed role."""
import hashlib
import json
from pathlib import Path
import pickle
import pytest


class UnreviewedContainer:
    pass


def test_installed_foundation_rejects_unreviewed_reconstruction(tmp_path):
    import torch
    from fairchem.core.units.mlip_unit.api.inference import MLIPInferenceCheckpoint
    from src.simulations.peregrine.foundation import load_foundation
    from peregrine.artifacts.store import canonical
    for i,config in enumerate([{'_target_':'os.system'}, {'head_cls':'unknown.Constructor'}, {'name':'${oc.env:UNREVIEWED}'}]):
        staged=tmp_path/('staged-'+str(i));staged.mkdir()
        saved=MLIPInferenceCheckpoint(model_config=config,model_state_dict={},ema_state_dict={},tasks_config={})
        checkpoint=staged/'checkpoint.pt';torch.save(saved,checkpoint)
        spec=dict(schema='anvil.uma-foundation/v1',sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),species=['Pt','O'],cutoff=6.0,task='oc20')
        raw=canonical(spec);identity=hashlib.sha256(raw).hexdigest();(staged/'foundation.json').write_bytes(raw)
        staged.rename(tmp_path/identity)
        with pytest.raises(ValueError,match='Unreviewed|interpolation'):load_foundation(tmp_path,identity,'cpu')


@pytest.mark.parametrize('payload,exception', [(UnreviewedContainer,pickle.UnpicklingError), (lambda:complex(1,2),ValueError)])
def test_installed_foundation_rejects_pickle_global_and_digest(tmp_path,payload,exception):
    import torch
    from src.simulations.peregrine.foundation import load_foundation
    from peregrine.artifacts.store import canonical
    checkpoint=tmp_path/'checkpoint.pt';torch.save(payload(),checkpoint)
    spec=dict(schema='anvil.uma-foundation/v1',sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),species=['Pt'],cutoff=6.0,task='oc20')
    raw=canonical(spec);identity=hashlib.sha256(raw).hexdigest();entry=tmp_path/identity;entry.mkdir();checkpoint.rename(entry/'checkpoint.pt');(entry/'foundation.json').write_bytes(raw)
    with pytest.raises(exception):load_foundation(tmp_path,identity,'cpu')
    (entry/'checkpoint.pt').write_bytes(b'damaged')
    with pytest.raises(ValueError,match='hash differs'):load_foundation(tmp_path,identity,'cpu')


def test_foundation_result_identity_comes_from_verified_spec(tmp_path, monkeypatch):
    import torch
    from fairchem.core.units.mlip_unit.api.inference import MLIPInferenceCheckpoint
    from src.simulations.peregrine.foundation import load_foundation
    from peregrine.artifacts.store import canonical
    import peregrine.train.factory as factory
    checkpoint = tmp_path / 'checkpoint.pt'
    torch.save(MLIPInferenceCheckpoint(model_config={}, model_state_dict={}, ema_state_dict={}, tasks_config={}), checkpoint)
    spec = dict(schema='anvil.uma-foundation/v1', sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(), species=['Pt'], cutoff=6.0, task='oc22')
    raw = canonical(spec)
    identity = hashlib.sha256(raw).hexdigest()
    entry = tmp_path / identity
    entry.mkdir()
    checkpoint.rename(entry / 'checkpoint.pt')
    (entry / 'foundation.json').write_bytes(raw)
    def build(model_spec, **kwargs):
        assert model_spec.backbone.kwargs['task_name'] == 'oc22'
        return torch.nn.Linear(1, 1)
    monkeypatch.setattr(factory, 'build_model_from_spec', build)
    _, metadata = load_foundation(tmp_path, identity, 'cpu')
    assert metadata['model_identity'] == dict(model=identity, checkpoint_sha256=spec['sha256'], task='oc22')
