"""Public Peregrine operations with caller-owned JSON evidence and artifacts."""
import json
from pathlib import Path

UNITS = dict(energy='eV', forces='eV/angstrom', stress='eV/angstrom^3', positions='angstrom', time='ps', temperature='K')


def validate_restart(request, root, atoms, saved):
    import numpy as np
    previous = json.loads((Path(root).parent / 'simulation.json').read_text())['request']
    for key in ('model', 'seed'):
        if previous[key] != request[key]:
            raise ValueError('Restart model/seed differs from request')
    if previous.get('system_ids') != request.get('system_ids') or [a.get('fixed_atoms', []) for a in previous['structures']] != [a.get('fixed_atoms', []) for a in request['structures']]:
        raise ValueError('Restart system order or fixed atoms differ from request')
    previous_identity = json.loads((Path(root).parent / 'simulation.json').read_text())['implementation']
    from src.simulations.identity import identity
    if previous_identity != identity("peregrine"):
        raise ValueError('Restart scientific implementation changed')
    for a, b in zip(atoms, saved, strict=True):
        if a.info != b.info or any(not np.array_equal(getattr(a, k), getattr(b, k)) for k in ('positions', 'numbers', 'cell', 'pbc')):
            raise ValueError('Restart geometry/metadata differs from request')


def prepare_neb(model, atoms, request, root):
    from peregrine_sim.simulation.neb import NEBSession
    options = {k:v for k,v in request['options'].items() if k != 'steps'}
    if request.get('restart'):
        session = NEBSession.load(Path(root).parent / 'checkpoint', model, model_id=request['model'], expected_digest=request['restart']['manifest_sha256'], expected_config=options)
        validate_restart(request, root, atoms, session.result.to_atoms())
        return session
    return NEBSession(model, atoms, model_id=request['model'], **options)


def extended(model, atoms, ids, request, output, session=None):
    from peregrine_sim import run_static
    from peregrine.artifacts.store import canonical
    operation, options = request['operation'], dict(request['options'])
    checkpoint = None
    if operation == 'neb':
        steps = options.pop('steps')
        result = session.run(steps, trajectory=output / 'neb-progress.xyz')
        checkpoint = dict(path='checkpoint', manifest_sha256=session.save(output / 'checkpoint'))
        physical = result.report()
        final = result.to_atoms()
        from ase.io import write
        write(output / 'images.xyz', final)
        report = dict(schema_version=1, engine='torch_sim', operation=operation, error=None, units=UNITS,
                      systems=[dict(system_id=i, steps=result.steps, converged=result.converged, termination=physical['termination']) for i in ids],
                      observables=[dict(system_id=i, energy=float(e), forces=f.tolist()) for i,e,f in zip(ids, result.energies, result.forces.split([len(a) for a in atoms]), strict=True)])
        (output / 'path-profile.json').write_bytes(canonical(physical))
    else:
        final = atoms
        if operation == 'frequencies':
            from peregrine_sim.analysis.neb import analyze_frequencies
            import numpy as np
            physical = []
            for i, a in enumerate(atoms):
                frequency = analyze_frequencies(model, a, **options)
                physical.append(frequency.report())
                np.savez(output / ('frequency-%d.npz' % i), **{key:getattr(frequency,key).cpu().numpy() for key in ('hessian', 'eigenvalues', 'modes', 'atom_indices')})
            (output / 'frequencies.json').write_bytes(canonical(physical))
        else:
            from peregrine_sim.sampling.generators import rattling
            from peregrine_sim.sampling.selectors import select_structures
            from peregrine_sim.io import state_to_atoms
            if len(atoms) != 1:
                raise ValueError('Sampling generation requires exactly one seed structure')
            n_select, strategy, descriptor = (options.pop(k) for k in ('n_select', 'strategy', 'descriptor'))
            seed = atoms[0].copy()
            mobile = [True] * len(seed)
            for constraint in seed.constraints:
                for index in constraint.get_indices():
                    mobile[index] = False
            seed.set_constraint()
            states = rattling(seed, mobile=mobile, seed=request['seed'], device=model.device, dtype=model.dtype, **options)
            candidates = [ids[0] + ':' + str(i) for i in range(states.n_systems)]
            if descriptor == 'backbone':
                from peregrine.descriptors import BackboneDescriptor
                descriptor = BackboneDescriptor(model.model, level='structure')
            _, selection = select_structures(states, n_select=n_select, model=model, strategy=strategy, descriptor=descriptor, identities=candidates, seed=request['seed'], return_result=True)
            physical = selection.to_dict()
            (output / 'selection.json').write_bytes(canonical(physical))
            all_atoms = state_to_atoms(states)
            for candidate in all_atoms:
                candidate.info.update(atoms[0].info)
                candidate.set_constraint(atoms[0].constraints)
                if any(not mobile[i] and (candidate.positions[i] != atoms[0].positions[i]).any() for i in range(len(seed))):
                    raise ValueError('Candidate generation moved fixed atoms')
            from ase.io import write
            write(output / 'candidates.xyz', all_atoms)
            final = [all_atoms[i] for i in selection.indices]
            ids = [candidates[i] for i in selection.indices]
        if final:
            report = run_static(model, final, system_ids=ids, return_result=True).report()
            report['operation'] = operation
        else:
            report = dict(schema_version=1, engine='torch_sim', operation=operation, error=None, units=UNITS, systems=[], observables=[])
    return final, report, checkpoint, physical, ids


