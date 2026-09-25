"""Independent periodic simulations through Peregrine; native checkpoints remain engine-owned."""
import json
import math
from pathlib import Path
import shutil
import sys


from src.simulations.peregrine.validation import validate_request_api
from src.simulations.peregrine.operations import extended, prepare_neb, validate_restart


def same_at_model_precision(expected, observed, dtype):
    """The fixed cell/constraints must equal the tensors actually passed to the model."""
    import numpy as np
    return np.array_equal(np.asarray(expected).astype(str(dtype).removeprefix('torch.')), observed)


def biased_model(model, request):
    """Construct only public provider objects; numerical history remains provider-owned."""
    from peregrine_sim.bias import BiasedModel, UmbrellaRestraint
    from peregrine_sim.colvars import Distance, Angle, Dihedral, StackedCV
    from peregrine_sim.simulation import BiasRecorder, metadynamics_bias
    spec = request["bias"]
    classes = {"distance": Distance, "angle": Angle, "dihedral": Dihedral}
    def build(c):
        if c['kind'] == 'position_projection':
            from peregrine_sim.colvars import PositionProjection
            return PositionProjection(c['atoms'][0], c['origin'], c['direction'])
        return classes[c['kind']](*c['atoms'])
    cv = StackedCV([build(c) for c in spec['cvs']])
    if spec["kind"] == "umbrella":
        bias = UmbrellaRestraint(cv, spec["centers"], spec["force_constants"])
    else:
        bias = metadynamics_bias(cv, height=spec["height"], sigma=spec["sigma"], pace=spec["pace"], temperature=request["options"]["temperature"], bias_factor=spec.get("bias_factor"))
    return BiasedModel(model, bias).eval(), BiasRecorder(system_ids=request["system_ids"])



def execute(request, registry_root, output):
    validate_request_api(request)
    import numpy as np
    import torch
    from ase import Atoms
    from peregrine.artifacts import ModelRegistry
    from peregrine.artifacts.store import canonical, file_digest
    from peregrine.interface.torchsim import PeregrineModel
    from peregrine_sim import relax, run_static
    from peregrine_sim.simulation import SimulationSession
    from torch_sim.io import state_to_atoms

    required = {"operation", "structures", "options", "seed", "model"}
    if not required <= set(request) or set(request) - required - {"system_ids", "restart", "bias", "device", "record_interval", "ensemble", "friction", "transport", "checkpoint_interval"}:
        raise ValueError("Unexpected worker request fields")
    operation, options = request["operation"], request["options"]
    restart = request.get("restart")
    if restart is not None and (operation not in {"md", "neb"} or set(restart) != {"manifest_sha256"}):
        raise ValueError("Unsupported restart request")
    from ase.constraints import FixAtoms
    atoms = []
    for row in request['structures']:
        fixed = row.get('fixed_atoms', [])
        if any(i >= len(row['numbers']) for i in fixed):
            raise ValueError('Fixed atom index outside structure')
        a = Atoms(**{k:v for k,v in row.items() if k != 'fixed_atoms'})
        if fixed:
            a.set_constraint(FixAtoms(indices=fixed))
        atoms.append(a)
    ids = request.get("system_ids", [str(i) for i in range(len(atoms))])
    if len(ids) != len(atoms) or any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError("Ordered unique system IDs required")
    if any(not len(a) or a.cell.rank != 3 or not np.isfinite(a.positions).all() or not np.isfinite(a.cell).all() or np.linalg.det(a.cell) == 0 for a in atoms):
        raise ValueError("This bridge requires finite structures with nonsingular cells")
    source = ModelRegistry(registry_root)
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    foundation = (Path(registry_root) / request['model'] / 'foundation.json').is_file()
    device = request.get('device', 'cpu')
    pending, seen = ([] if foundation else [request["model"]]), set()
    while pending:
        version = pending.pop()
        if version in seen:
            continue
        manifest = source.store.read(version)
        meta = manifest["metadata"]
        if meta.get("foundation", {}).get("files"):
            raise ValueError("External foundation files require a deployment adapter")
        if meta["kind"] == "deployment":
            architecture = meta["architecture"]
            recipe = architecture.get("recipe", {})
            if architecture.get("kind") != "hybrid" or recipe.get("backbone", {}).get("module") != "peregrine.backbones.native" or recipe.get("corrections") or recipe.get("projection"):
                raise ValueError("This CPU bridge accepts native hybrid deployments without corrections/projection only")
        seen.add(version)
        pending.extend(meta[k] for k in ("parent", "dataset") if meta.get(k))
    if foundation:
        from src.simulations.peregrine.foundation import load_foundation
        potential, meta = load_foundation(registry_root, request['model'], device)
    else:
        meta = source.manifest(request['model'])
    model_identity = meta.get('model_identity', dict(model=request['model']))
    for a in atoms:
        if set(a.info) - {"charge", "spin", "head", "fidelity_idx"}:
            raise ValueError("Unsupported structure metadata; refusing silent loss")
        for key, permitted in meta["target"]["selectors"].items():
            value = a.info.get(key, 0)
            if type(value) is not int or value not in permitted:
                raise ValueError("Prediction selector outside deployment target")
        for key in ("charge", "spin"):
            if key in a.info and (type(a.info[key]) not in (int, float) or not np.isfinite(a.info[key])):
                raise ValueError("Charge and spin must be finite scalars")
    if any(set(a.info) != set(atoms[0].info) for a in atoms):
        raise ValueError("Batch structures must declare the same metadata keys")
    # Service reconstruction validates every fitted quantifier, including sidecars.
    if not foundation:
        source.service(request["model"])
        potential = source.load(request["model"], device=device)
    if any(set(a.get_chemical_symbols()) - set(potential.species) for a in atoms):
        raise ValueError("Structure species outside model vocabulary")
    model = PeregrineModel(potential, compute_forces=True, compute_stress=not foundation, device=device)
    recorder = None
    if "bias" in request:
        model, recorder = biased_model(model, request)
    torch.manual_seed(request["seed"])
    checkpoint, prior_steps, session = None, 0, None
    if operation == 'neb':
        session = prepare_neb(model, atoms, request, registry_root)
    if operation == "md":
        from torch_sim.units import UnitSystem
        config = dict(mode=request.get("ensemble", "nve"), temperature=float(options["temperature"]), timestep=float(options["timestep"]), force_tol=None, seed=request["seed"], options={"gamma": request["friction"] / UnitSystem.metal.time} if "friction" in request else {}, init_options={})
        if restart is None:
            initial = atoms
            session = SimulationSession(model, initial, model_id=request["model"], system_ids=ids, **config)
        else:
            digest = restart["manifest_sha256"]
            if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("An externally recorded checkpoint digest is required")
            session = SimulationSession.load(Path(registry_root).parent / "checkpoint", model, model_id=request["model"], expected_digest=digest, expected_config=dict(config, engine="torch_sim"))
            saved = state_to_atoms(session.final, system_extras_map={k: k for k in atoms[0].info})
            saved = saved if isinstance(saved, list) else [saved]
            if [s["system_id"] for s in session.result().report()["systems"]] != ids:
                raise ValueError("Restart system order differs from request")
            prior_steps = session.result().report()["systems"][0]["steps"]
            if any(s["steps"] != prior_steps for s in session.result().report()["systems"]):
                raise ValueError("Restart system step counts differ")
            validate_restart(request, registry_root, atoms, saved)
    output.mkdir()  # Identity and restart validation precede output side effects.
    staged = ModelRegistry(output / "deployment")
    if foundation:
        from src.simulations.peregrine.foundation import copy_foundation
        copy_foundation(registry_root, request['model'], output / 'deployment')
    for version in sorted(seen):
        shutil.copytree(source.store.path(version), staged.store.path(version))
    frames, physical = dict(records=[], count=0), None
    def observe(state, step):
        if recorder:
            recorder.observe(state, step)
        if step % request.get('record_interval', 1) == 0:
            from ase.io import write
            geometries = state_to_atoms(state, system_extras_map={k:k for k in atoms[0].info})
            for identity, initial, geometry in zip(ids, atoms, geometries, strict=True):
                geometry.calc = None
                geometry.info.update(physical_step=step, system_id=identity, label_semantics='geometry_only')
                geometry.set_constraint(initial.constraints)
            write(output / 'trajectory.xyz', geometries, append=True)
            temperatures = (state.calc_kT()/UnitSystem.metal.temperature).detach().cpu().tolist()
            energies = state.energy.detach().cpu().tolist()
            with (output/'thermodynamics.jsonl').open('a') as stream:
                stream.write(json.dumps(dict(step=step, time_ps=step*float(options['timestep']), temperature_K=temperatures, energy_eV=energies), allow_nan=False)+'\n')
            if frames['count'] < 1000:
                for i, identity in enumerate(ids):
                    energy = recorder.records[-1]['base_energy'][i] if recorder else energies[i]
                    frames['records'].append(dict(ordinal=frames['count'], coordinate=step * float(options['timestep']), energy=energy, data=dict(step=step, frame_offset=frames['count'] * len(ids) + i, system_id=identity, temperature_K=temperatures[i], geometry_file='trajectory.xyz')))
            frames['count'] += 1
    if operation in {'neb', 'frequencies', 'sampling'}:
        final, report, checkpoint, physical, ids = extended(model, atoms, ids, request, output, session)
        if operation == 'sampling':
            atoms = final
    elif operation == "md":
        if restart is None and recorder is None:observe(session.final, 0)
        remaining = options['steps']
        interval = request.get('checkpoint_interval', remaining)
        while remaining:
            chunk = min(interval, remaining)
            result = session.run(chunk, observer=observe)
            remaining -= chunk
            if remaining:
                session.save(output / ('recovery-%d' % result.report()['systems'][0]['steps']))
        checkpoint = dict(path="checkpoint", manifest_sha256=session.save(output / "checkpoint"))
    else:
        result = (run_static if operation == "static" else relax)(model, atoms, system_ids=ids, return_result=True, **options)
    if operation not in {"neb", "frequencies", "sampling"}:
        report = result.report()
    if operation == 'md':
        measured = (result.final.calc_kT() / UnitSystem.metal.temperature).detach().cpu().tolist()
        if len(measured) != len(ids) or not np.isfinite(measured).all():
            raise ValueError('Invalid measured kinetic temperature')
        physical = dict(ensemble=config['mode'], target_temperature_K=options['temperature'], measured_temperature_K=measured, temperature_basis='kinetic energy with native constrained degrees of freedom')
        if 'transport' in request:
            from src.utils.analysis.transport import analyze
            physical['transport'] = analyze(output/'trajectory.xyz', request['transport'], timestep_ps=options['timestep'], temperature_K=options['temperature'], output_dir=output)
            (output/'transport.json').write_bytes(canonical(physical['transport']))
        report['thermal'] = physical
    if recorder is not None:
        if [row["step"] for row in recorder.records] != list(range(prior_steps + 1, prior_steps + options["steps"] + 1)):
            raise ValueError("Missing or duplicated completed-step bias observations")
        for row in recorder.records:
            for key, shape in {"base_energy": (len(ids),), "bias_energy": (len(ids),), "total_energy": (len(ids),), "cv_0": (len(ids), len(request["bias"]["cvs"]))}.items():
                value = np.asarray(row[key], dtype=float)
                if value.shape != shape or not np.isfinite(value).all():
                    raise ValueError("Invalid bias observation shape or values")
            if row["system_ids"] != ids or row["energy_unit"] != "eV":
                raise ValueError("Bias observation identity/units mismatch")
        (output / "bias-observations.json").write_bytes(canonical(dict(schema="atomistic-skills.bias-observations/v1", configuration=request["bias"], records=recorder.records)))
    if report["units"] != dict(energy="eV", forces="eV/angstrom", stress="eV/angstrom^3", positions="angstrom", time="ps", temperature="K"):
        raise ValueError("Unqualified simulation units")
    if report["schema_version"] != 1 or report["engine"] != "torch_sim" or report["operation"] != operation or report["error"] is not None or [s["system_id"] for s in report["systems"]] != ids or [v["system_id"] for v in report["observables"]] != ids:
        raise ValueError("Unsupported, failed or unordered simulation report")
    if operation not in {"neb", "frequencies", "sampling"}:
        final = atoms if operation == "static" else state_to_atoms(result.final, system_extras_map={k: k for k in atoms[0].info})
    final = final if isinstance(final, list) else [final]
    rows = []
    for i, (before, after) in enumerate(zip(atoms, final, strict=True)):
        if before.info != after.info or not np.array_equal(before.numbers, after.numbers) or not np.array_equal(before.pbc, after.pbc):
            raise ValueError("Worker changed metadata, atom order/species or PBC")
        if not np.isfinite(after.positions).all() or not same_at_model_precision(before.cell, after.cell, model.dtype):
            raise ValueError("Invalid final geometry or unexpected cell change")
        fixed = before.constraints[0].get_indices() if before.constraints else []
        if len(fixed) and not same_at_model_precision(before.positions[fixed], after.positions[fixed], model.dtype):
            raise ValueError('Fixed atoms moved')
        geometry = dict(numbers=after.numbers.tolist(), positions=after.positions.tolist(), cell=after.cell.array.tolist(), pbc=after.pbc.tolist(), info=before.info)
        if len(fixed):
            geometry["fixed_atoms"] = list(map(int, fixed))
        values, outcome = report["observables"][i], report["systems"][i]
        if type(outcome["steps"]) is not int or outcome["steps"] < 0:
            raise ValueError("Invalid actual step count")
        if operation in {"relax", "neb"}:
            if type(outcome["converged"]) is not bool or (operation == "relax" and outcome["steps"] > options["max_steps"]) or outcome["termination"] != ("converged" if outcome["converged"] else "step_limit"):
                raise ValueError("Invalid relaxation outcome")
        elif outcome["converged"] is not None or outcome["termination"] != "completed" or (operation == "static" and outcome["steps"] != 0) or (operation == "md" and outcome["steps"] != prior_steps + options["steps"]):
            raise ValueError("Invalid static/MD outcome")
        for key, shape in {"energy": (), "forces": (len(after), 3), "stress": (3, 3)}.items():
            if key == "stress" and key not in values and (operation != "static" or foundation):
                continue  # Native motion reports need not retain stress.
            value = np.asarray(values.get(key), dtype=float)
            if value.shape != shape or not np.isfinite(value).all():
                raise ValueError(f"Missing, null, nonfinite or incorrectly shaped {key}: {value.shape}")
        data = dict(values, geometry={k:v for k,v in geometry.items() if k in ("info", "fixed_atoms")}, units=report["units"], stress_sign="tensile-positive",
                    model_identity=model_identity,
                    model=request["model"], operation=operation, options=options, seed=request["seed"], label_kind="ml_prediction", convergence=outcome["converged"], outcome=outcome, engine="torch_sim", checkpoint=checkpoint)
        if operation == 'md':
            data.update(ensemble=config['mode'], target_temperature_K=options['temperature'], measured_temperature_K=measured[i], transport=physical.get('transport'))
        if recorder is not None:
            observation = recorder.records[-1]
            data.update(label_kind="biased_trajectory", energy_basis="biased_total", forces_basis="biased_total", base_energy=observation['base_energy'][i], bias_energy=observation['bias_energy'][i], total_energy=observation['total_energy'][i], bias_configuration=request["bias"], bias_observations="bias-observations.json")
        after.set_constraint(before.constraints)
        # The immutable worker export remains independently readable; DB CalcNode has references only.
        rows.append(dict(data, geometry=geometry))
    from src.simulations.identity import identity
    document = dict(schema="atomistic-skills.simulation/v1", request=request, report=report,
        checkpoint=checkpoint, results=rows, dependencies=sorted(seen), implementation=identity("peregrine"),
        physical=physical, frames=frames, model_identity=model_identity)
    (output / "simulation.json").write_bytes(canonical(document))
    return document
