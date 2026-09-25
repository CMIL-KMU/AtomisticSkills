"""Reuse upstream diffusion physics while checking units and provisional reporting."""

import json
import numpy as np
import pytest

pytest.importorskip("pymatgen.analysis.diffusion.analyzer")
from ase import Atoms
from ase.io import write
from src.utils.analysis.transport import analyze


def trajectory(tmp_path, moving=True):
    random = np.random.default_rng(9)
    coordinates = np.zeros((33, 3))
    coordinates[:, 0] = np.arange(33)
    frames = []
    for step in range(801):
        if moving:
            coordinates[:32] += random.normal(0, 0.08, (32, 3))
        a = Atoms("Li32Cl", positions=coordinates, cell=[60, 60, 60], pbc=True)
        a.info.update(physical_step=step * 10, system_id="sample")
        frames.append(a)
    path = tmp_path / "trajectory.xyz"
    write(path, frames)
    return path


def test_ne_units_and_no_claim_of_convergence(tmp_path):
    from scipy.constants import elementary_charge, Boltzmann

    path = trajectory(tmp_path)
    result = analyze(
        path,
        dict(species="Li", charge=1, equilibration_ps=1, min_observations=30),
        timestep_ps=0.001,
        temperature_K=600,
        output_dir=tmp_path,
    )
    assert result["status"] == "provisional" and result["production_duration_ps"] == 7
    expected = (
        result["diffusivity_cm2_s"]
        * 1e-4
        * (32 / (60**3 * 1e-30))
        * elementary_charge**2
        / (Boltzmann * 600)
        * 10
    )
    assert result["conductivity_NE_mS_cm"] == pytest.approx(expected)
    assert (tmp_path / "msd.png").is_file() and (tmp_path / "msd.csv").is_file()
    json.dumps(result, allow_nan=False)


def test_stationary_atoms_do_not_become_positive_diffusion(tmp_path):
    result = analyze(
        trajectory(tmp_path, False),
        dict(species="Li", charge=1, equilibration_ps=0, min_observations=30),
        timestep_ps=0.001,
        temperature_K=600,
    )
    assert (
        result["diffusivity_cm2_s"] is None and result["conductivity_NE_mS_cm"] is None
    )


@pytest.mark.parametrize("smoothed", [False, "max"])
def test_pymatgen_parity_and_explicit_fit_window(tmp_path, smoothed):
    from ase.io import read
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer

    path = trajectory(tmp_path)
    cfg = dict(species="Li", charge=1, equilibration_ps=0, smoothed=smoothed)
    result = analyze(path, cfg, timestep_ps=0.001, temperature_K=600)
    structures = [AseAtomsAdaptor.get_structure(a) for a in read(path, ":")]
    reference = DiffusionAnalyzer.from_structures(
        structures, "Li", 600, time_step=1, step_skip=10, smoothed=smoothed, min_obs=30
    )
    assert result["diffusivity_cm2_s"] == pytest.approx(
        reference.diffusivity, rel=1e-10
    )
    assert result["fit_standard_error_cm2_s"] == pytest.approx(
        reference.diffusivity_std_dev, rel=1e-10
    )
    window = analyze(
        path,
        dict(cfg, fit_start_ps=2, fit_end_ps=4),
        timestep_ps=0.001,
        temperature_K=600,
    )
    assert window["fit_lag_ps"][0] >= 2 and window["fit_lag_ps"][1] <= 4


def test_non_lithium_charge_scaling_and_explicit_generic_timing(tmp_path):
    from ase.io import read

    path = trajectory(tmp_path)
    frames = read(path, ":")
    for a in frames:
        a.numbers[a.numbers == 3] = 11
        a.info.clear()
    write(path, frames)
    cfg = dict(species="Na", charge=1, equilibration_ps=0)
    with pytest.raises(ValueError, match="frame_interval"):
        analyze(path, cfg, temperature_K=600)
    first = analyze(path, cfg, temperature_K=600, frame_interval_fs=10)
    second = analyze(path, dict(cfg, charge=2), temperature_K=600, frame_interval_fs=10)
    assert first["species"] == "Na"
    assert second["conductivity_NE_mS_cm"] == pytest.approx(
        4 * first["conductivity_NE_mS_cm"]
    )


def test_negative_raw_slope_is_not_positive_floor(tmp_path):
    from ase.io import read

    path = trajectory(tmp_path, False)
    frames = read(path, ":")
    # Ions move outward then back; the selected return interval has negative slope.
    for i, a in enumerate(frames):
        a.positions[:32, 1] = min(i, 800 - i) * 0.002
    write(path, frames)
    result = analyze(
        path,
        dict(species="Li", charge=1, equilibration_ps=0, fit_start_ps=5),
        temperature_K=600,
        timestep_ps=0.001,
    )
    assert result["fit_slope_angstrom2_fs"] < 0
    assert result["diffusivity_cm2_s"] is None


def test_changed_cell_and_invalid_configuration_fail(tmp_path):
    from ase.io import read
    from src.utils.analysis.transport import configuration

    for change in [
        dict(smoothed=True),
        dict(charge=0),
        dict(fit_start_ps=5, fit_end_ps=2),
        dict(equilibration_ps=float("nan")),
    ]:
        with pytest.raises(ValueError):
            configuration(
                dict({"species": "Li", "charge": 1, "equilibration_ps": 0}, **change)
            )
    path = trajectory(tmp_path)
    frames = read(path, ":")
    frames[-1].cell[0, 0] += 1
    write(path, frames)
    with pytest.raises(ValueError, match="fixed cell"):
        analyze(
            path,
            dict(species="Li", charge=1, equilibration_ps=0),
            temperature_K=600,
            timestep_ps=0.001,
        )


def test_conflicting_timing_and_max_lag_resolution_are_rejected(tmp_path):
    from ase.io import read

    path = trajectory(tmp_path, False)
    frames = read(path, ":")
    for frame in frames:
        frame.info.clear()
    write(path, frames)
    cfg = dict(species="Li", charge=1, equilibration_ps=0)
    with pytest.raises(ValueError, match="not both"):
        analyze(path, cfg, temperature_K=600, timestep_ps=0.001, frame_interval_fs=10)
    with pytest.raises(ValueError, match="spacing <= 1 ps"):
        analyze(
            path, dict(cfg, smoothed="max"), temperature_K=600, frame_interval_fs=2000
        )
