"""Data leakage, canonical labels and metric units at the Peregrine skill boundary."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from pymatgen.io.ase import AseAtomsAdaptor

SCRIPT = Path(__file__).resolve().parents[1] / '.agents/skills/ml-mace-finetune/scripts/peregrine_finetune.py'
spec = importlib.util.spec_from_file_location('peregrine_finetune', SCRIPT)
script = importlib.util.module_from_spec(spec)
spec.loader.exec_module(script)


def dataset(tmp_path):
    rows = []
    for i, split in enumerate(['train', 'validation', 'test']):
        atoms = bulk('Pt', 'fcc', a=3.9 + i * .01, cubic=True)
        rows.append(dict(id=str(i), split=split, structure=AseAtomsAdaptor.get_structure(atoms).as_dict(),
                         energy=-20., forces=np.zeros((4,3)).tolist(), stress=np.zeros((3,3)).tolist()))
    doc = dict(units=dict(energy='eV', forces='eV/angstrom', stress='eV/angstrom3', stress_sign='tensile'), structures=rows)
    path = tmp_path / 'data.json'
    path.write_text(json.dumps(doc))
    return path, doc


@pytest.mark.parametrize('bad', ['units', 'duplicate', 'nonfinite', 'shape', 'split', 'stress'])
def test_reject_ambiguous_or_leaking_labels(tmp_path, bad):
    path, doc = dataset(tmp_path)
    if bad == 'units': doc['units']['stress'] = 'kbar'
    if bad == 'duplicate': doc['structures'][2]['structure'] = doc['structures'][0]['structure']
    if bad == 'nonfinite': doc['structures'][0]['energy'] = float('nan')
    if bad == 'shape': doc['structures'][0]['forces'] = [[0., 0., 0.]]
    if bad == 'split': doc['structures'][2]['split'] = 'train'
    if bad == 'stress': doc['structures'][0]['stress'][0][1] = 1.
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError): script.load_dataset(path)


def test_explicit_split_matches_native_datamodule(tmp_path):
    torch = pytest.importorskip('torch')
    pytest.importorskip('peregrine')
    from peregrine.train.datamodule import PeregrineDataModule
    from peregrine.train.unified_config import UnifiedTrainingConfig
    from ase.io import write
    path, _ = dataset(tmp_path)
    rows, frames = script.load_dataset(path)
    order, fraction = script.training_order(rows, 57)
    assert set(order) == {0, 1}  # Test frame never reaches the trainer.
    source = tmp_path / 'training.extxyz'
    write(source, [frames[i] for i in order])
    config = UnifiedTrainingConfig.from_dict(dict(model=dict(strategy='native_passthrough',backbone=dict(type='mace', species=['Pt'])),
        data=dict(sources=[dict(path=str(source))], targets=['energy','forces','stress'],
                  split=dict(val_fraction=fraction,seed=57))))
    dm = PeregrineDataModule(config.data,species=['Pt'],cutoff=5.,dtype=torch.float32)
    dm.setup()
    assert len(dm._train) == len(dm._val) == 1
    assert np.isclose(float(dm._train[0]['cell'].reshape(3,3)[0,0]),frames[0].cell[0,0])
    assert np.isclose(float(dm._val[0]['cell'].reshape(3,3)[0,0]),frames[1].cell[0,0])


def test_metric_normalization(tmp_path):
    from ase.units import GPa
    path, _ = dataset(tmp_path)
    rows, frames = script.load_dataset(path)
    predictions = [dict(energy=r['energy']+.004,forces=np.ones((4,3)).tolist(),
                        stress=(np.ones((3,3))*GPa).tolist()) for r in rows]
    for metrics in script.metrics(rows,frames,predictions).values():
        assert metrics['energy_meV_per_atom']['mae'] == pytest.approx(1.)
        assert metrics['forces_eV_per_angstrom']['mae'] == pytest.approx(1.)
        assert metrics['stress_GPa']['mae'] == pytest.approx(1.)
