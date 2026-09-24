#!/usr/bin/env python
# Env: peregrine-mace
"""Fine-tune a trusted local MACE checkpoint using Peregrine and explicit E/F/stress splits."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from ase.stress import full_3x3_to_voigt_6_stress
from ase.units import GPa
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from src.utils.analysis.metrics import evaluate_metrics


def checksum(path):
    """Bind exact dataset/checkpoint bytes."""
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def provider_identity():
    """Record installed provider bytes even when development builds share a version."""
    import peregrine
    root = Path(peregrine.__file__).parent
    files = {p.relative_to(root).as_posix(): checksum(p) for p in sorted(root.rglob('*.py'))}
    return dict(files=files, sha256=hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest())


def load_dataset(path):
    """Reject missing labels, duplicate geometries and ambiguous unit conventions."""
    doc = json.loads(Path(path).read_text())
    if doc['units'] != dict(energy='eV', forces='eV/angstrom', stress='eV/angstrom3', stress_sign='tensile'):
        raise ValueError('Explicit canonical E/F/stress units required; convert raw labels once upstream')
    rows, frames, ids, geometries = doc['structures'], [], set(), set()
    for row in rows:
        if not isinstance(row['id'], str) or not row['id'] or row['id'] in ids:
            raise ValueError('Unique nonempty structure IDs required')
        if row['split'] not in {'train', 'validation', 'test'}:
            raise ValueError('Explicit train/validation/test assignment required')
        atoms = AseAtomsAdaptor.get_atoms(Structure.from_dict(row['structure']))
        if (not len(atoms) or not all(atoms.pbc) or atoms.get_volume() <= 0
                or not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell.array).all()):
            raise ValueError('Nonempty periodic cells required for stress training')
        geometry = json.dumps(dict(z=atoms.numbers.tolist(), p=np.round(atoms.positions, 10).tolist(),
                                   c=np.round(atoms.cell.array, 10).tolist()), sort_keys=True)
        if geometry in geometries:
            raise ValueError('Duplicate geometry would contaminate the split')
        e, f, s = np.asarray(row['energy']), np.asarray(row['forces']), np.asarray(row['stress'])
        if e.shape != () or f.shape != (len(atoms), 3) or s.shape != (3, 3):
            raise ValueError('Energy/force/stress shape mismatch')
        if not all(np.isfinite(a).all() for a in (e, f, s)) or not np.allclose(s, s.T, atol=1e-10):
            raise ValueError('Finite energy/forces and symmetric stress required')
        atoms.calc = SinglePointCalculator(atoms, energy=float(e), forces=f, stress=s)
        atoms.info['structure_id'] = row['id']
        frames.append(atoms); ids.add(row['id']); geometries.add(geometry)
    if any(not any(row['split'] == split for row in rows) for split in ('train', 'validation', 'test')):
        raise ValueError('All three disjoint splits must be nonempty')
    return rows, frames


def training_order(rows, seed):
    """Reorder explicit splits to match Peregrine's native seeded DataModule split."""
    import torch
    train = [i for i, row in enumerate(rows) if row['split'] == 'train']
    valid = [i for i, row in enumerate(rows) if row['split'] == 'validation']
    n = len(train) + len(valid)
    positions = set(torch.randperm(n, generator=torch.Generator().manual_seed(seed)).tolist()[:len(valid)])
    train_iter, valid_iter = iter(train), iter(valid)
    return [next(valid_iter if i in positions else train_iter) for i in range(n)], len(valid) / n


def predict(model, frames, device):
    """Use Peregrine's existing ASE calculator for conserving E/F/stress."""
    import torch
    from peregrine.interface.ase import PeregrineCalculator
    calculator = PeregrineCalculator(model, device=device, dtype=torch.float32)
    predictions = []
    for frame in frames:
        atoms = frame.copy(); atoms.calc = calculator
        values = dict(energy=float(atoms.get_potential_energy()), forces=atoms.get_forces().tolist(),
                      stress=atoms.get_stress(voigt=False).tolist())
        if not all(np.isfinite(np.asarray(value)).all() for value in values.values()):
            raise ValueError('Model produced nonfinite predictions')
        predictions.append(values)
    return predictions


def metrics(rows, frames, predictions):
    """Reuse the public benchmark metrics with explicit component/atom normalization."""
    result = {}
    for split in ('train', 'validation', 'test'):
        indices = [i for i, row in enumerate(rows) if row['split'] == split]
        energy = evaluate_metrics([predictions[i]['energy'] / len(frames[i]) * 1000 for i in indices],
                                  [rows[i]['energy'] / len(frames[i]) * 1000 for i in indices])
        forces = evaluate_metrics(np.concatenate([predictions[i]['forces'] for i in indices]),
                                  np.concatenate([rows[i]['forces'] for i in indices]))
        stress = evaluate_metrics([full_3x3_to_voigt_6_stress(np.asarray(predictions[i]['stress'])) / GPa for i in indices],
                                  [full_3x3_to_voigt_6_stress(np.asarray(rows[i]['stress'])) / GPa for i in indices])
        result[split] = dict(structures=len(indices), energy_meV_per_atom=energy,
                             forces_eV_per_angstrom=forces, stress_GPa=stress)
    return result


def run(args):
    import torch
    from lightning.pytorch.loggers import CSVLogger
    from peregrine.train.driver import train_model
    from peregrine.train.factory import build_base_model
    from peregrine.train.trainer import PotentialTrainer
    from peregrine.train.unified_config import UnifiedTrainingConfig

    provider = provider_identity()
    if args.provider_sha256 and provider['sha256'] != args.provider_sha256:
        raise ValueError('Installed Peregrine source differs from the admitted provider')
    rows, frames = load_dataset(args.dataset)
    foundation = args.foundation.resolve(strict=True)
    if not 1 <= args.epochs <= 10000 or not 0 < args.lr <= 0.1 or args.batch_size < 1 or args.seed < 0:
        raise ValueError('Invalid training controls')
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('Use an empty output directory')
    before_hash, data_hash = checksum(foundation), checksum(args.dataset)
    order, fraction = training_order(rows, args.seed)
    write(output / 'train-validation.extxyz', [frames[i] for i in order])
    kwargs = dict(dtype='float32')
    if args.head: kwargs['head'] = args.head
    config = dict(model=dict(strategy='native_passthrough', backbone=dict(type='mace', source='pretrained',
        name_or_path=str(foundation), species=sorted(set(symbol for f in frames for symbol in f.get_chemical_symbols())),
        kwargs=kwargs), scale_shift=dict(source='none')),
        data=dict(sources=[dict(path=str(output / 'train-validation.extxyz'))], targets=['energy', 'forces', 'stress'],
                  split=dict(val_fraction=fraction, seed=args.seed), batch_size=args.batch_size),
        training=dict(seed=args.seed, max_epochs=args.epochs, optimizer=dict(type='adam', lr=args.lr),
            scheduler=dict(type='linear_warmup_cosine', kwargs=dict(warmup_epochs=min(10,args.epochs-1), cycle_epochs=args.epochs, num_cycles=1)),
            losses=[dict(type='mse',key='energy',weight=1.0),dict(type='mse',key='forces',weight=10.0),
                    dict(type='mse',key='stress',weight=1.0,scale=1/GPa)],
            metrics=[dict(type='mae',key=k) for k in ('energy','forces','stress')],
            save_top_k=1, early_stopping_patience=20, use_ema=False, log_every_n_steps=1))
    if args.trainable == 'readout':
        config['model']['finetune'] = [dict(name='selective',options=dict(trainable_modules=[r'\.readouts\.']))]
    cfg = UnifiedTrainingConfig.from_dict(config)
    save(output / 'training-config.json', config)
    model = build_base_model(cfg.model, dtype=torch.float32)
    trainable = {n:p.numel() for n,p in model.named_parameters() if p.requires_grad}
    if not trainable: raise ValueError('No trainable parameters')
    initial = {n:p.detach().cpu().clone() for n,p in model.named_parameters() if p.requires_grad}
    before = predict(model, frames, args.device)
    save(output / 'baseline.json', dict(predictions=before, metrics=metrics(rows,frames,before)))
    del model
    if args.device == 'cuda': torch.cuda.empty_cache()
    model, trainer = train_model(cfg, dtype=torch.float32, accelerator='gpu' if args.device=='cuda' else 'cpu',
        devices=1, default_root_dir=str(output / 'training'), logger=CSVLogger(str(output / 'training'),name='history'),
        enable_progress_bar=False, enable_model_summary=False, inference_mode=False)
    best = trainer.checkpoint_callback.best_model_path
    if not best: raise ValueError('Validation did not select a checkpoint')
    selected = PotentialTrainer.load_from_checkpoint(best, map_location='cpu', model=model).model
    changed = [n for n,p in selected.named_parameters() if n in initial and not torch.equal(initial[n],p.detach().cpu())]
    if not changed: raise ValueError('Training did not change any selected parameter')
    after = predict(selected, frames, args.device)
    np.savez(output / 'weights.npz', **{n:p.detach().cpu().numpy() for n,p in selected.state_dict().items()})
    selected.export_native(str(output / 'finetuned.model'))
    # Reload the portable export through the same public loader, independently of
    # Lightning's in-memory model/checkpoint, and compare all three observables.
    config['model']['backbone']['name_or_path'] = str(output / 'finetuned.model')
    restored = build_base_model(UnifiedTrainingConfig.from_dict(config).model, dtype=torch.float32)
    replay = predict(restored, frames, args.device)
    roundtrip = {key: max(float(np.max(np.abs(np.asarray(a[key])-np.asarray(b[key]))))
                          for a,b in zip(after,replay)) for key in ('energy','forces','stress')}
    # Float32 GPU derivative reductions can differ near zero after a reload.
    # These absolute physical tolerances are recorded, rather than relying on
    # a relative tolerance that becomes undefined for a zero force component.
    tolerances = dict(energy=1e-3, forces=1e-4, stress=1e-6)
    save(output / 'export-validation.json', dict(max_abs=roundtrip, tolerances=tolerances))
    for key in roundtrip:
        if roundtrip[key] > tolerances[key]:
            raise ValueError('Native model export changed '+key)
    if checksum(foundation) != before_hash or checksum(args.dataset) != data_hash:
        raise ValueError('Input changed during training')
    save(output / 'predictions.json', [dict(id=row['id'],split=row['split'],target={k:row[k] for k in ('energy','forces','stress')},
        before=before[i],after=after[i]) for i,row in enumerate(rows)])
    report = dict(schema='peregrine-mace-finetune/v1',dataset_sha256=data_hash,foundation_sha256=before_hash,
        checkpoint_sha256=checksum(output/'weights.npz'),native_checkpoint_sha256=checksum(output/'finetuned.model'),
        seed=args.seed,training_order=[rows[i]['id'] for i in order],test_ids=[r['id'] for r in rows if r['split']=='test'],
        before=metrics(rows,frames,before),after=metrics(rows,frames,after),trainable_parameters=sum(trainable.values()),
        changed_parameter_tensors=len(changed),best_validation_loss=float(trainer.checkpoint_callback.best_model_score),
        completed_epochs=trainer.current_epoch,export_roundtrip_max_abs=roundtrip,export_roundtrip_tolerances=tolerances,
        provider=provider,units=dict(energy='eV',forces='eV/angstrom',stress='eV/angstrom3',stress_sign='tensile'),
        versions={n:importlib.metadata.version(n) for n in ('torch','mace-torch','peregrine-pot','ase','numpy')})
    save(output / 'report.json', report)
    print(json.dumps(report,allow_nan=False))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset',type=Path,required=True)
    parser.add_argument('--foundation',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--seed',type=int,required=True)
    parser.add_argument('--epochs',type=int,required=True)
    parser.add_argument('--lr',type=float,required=True)
    parser.add_argument('--batch-size',type=int,required=True)
    parser.add_argument('--trainable',choices=['readout','all'],default='readout')
    parser.add_argument('--device',choices=['cuda','cpu'],default='cuda')
    parser.add_argument('--head')
    parser.add_argument('--provider-sha256',help='Optional expected installed Peregrine source SHA256')
    run(parser.parse_args())


if __name__=='__main__': main()
