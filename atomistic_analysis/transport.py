"""Fixed-cell ionic self diffusion using pymatgen; explicit timing and fit policy."""

import hashlib
import math
from importlib.metadata import version
from pathlib import Path


def configuration(config: dict) -> dict:
    """Resolve and validate defaults without guessing ion charge or trajectory timing."""
    allowed = {
        "species",
        "charge",
        "equilibration_ps",
        "min_observations",
        "smoothed",
        "fit_start_ps",
        "fit_end_ps",
    }
    if set(config) - allowed or not {"species", "charge", "equilibration_ps"} <= set(
        config
    ):
        raise ValueError(
            "Explicit species, charge and equilibration_ps required; unknown settings rejected"
        )
    cfg = dict(min_observations=30, smoothed=False, fit_start_ps=None, fit_end_ps=None)
    cfg.update(config)
    if not isinstance(cfg["species"], str) or not cfg["species"]:
        raise ValueError("species must be an element symbol")
    if (not isinstance(cfg["charge"], int) or isinstance(cfg["charge"], bool)) or cfg[
        "charge"
    ] == 0:
        raise ValueError("charge must be a nonzero integer in elementary-charge units")
    if cfg["smoothed"] is not False and cfg["smoothed"] != "max":
        raise ValueError('smoothed must be false or "max"')
    if (
        not isinstance(cfg["min_observations"], int)
        or isinstance(cfg["min_observations"], bool)
    ) or cfg["min_observations"] < 3:
        raise ValueError("min_observations must be an integer >= 3")
    for key in ("equilibration_ps", "fit_start_ps", "fit_end_ps"):
        value = cfg[key]
        if value is None and key != "equilibration_ps":
            continue
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    if cfg["fit_end_ps"] is not None and cfg["fit_end_ps"] <= (
        cfg["fit_start_ps"] or 0
    ):
        raise ValueError("fit_end_ps must exceed fit_start_ps")
    return cfg


def analyze(
    path,
    config: dict,
    *,
    temperature_K: float,
    timestep_ps=None,
    frame_interval_fs=None,
    output_dir=None,
) -> dict:
    """Analyze ASE-readable fixed-cell frames; write plots only to an explicit directory.

    Physical-step metadata requires timestep_ps (integrator step). Generic trajectories
    require frame_interval_fs (saved-frame spacing). Never infer timing from filenames.
    Results remain provisional. Overlapping MSD origins are not independent samples.
    """
    import numpy as np
    from ase.io import read
    from ase.data import atomic_numbers
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.analysis.diffusion.analyzer import (
        DiffusionAnalyzer,
        get_diffusivity_from_msd,
    )
    from scipy.constants import elementary_charge, Boltzmann

    cfg = configuration(config)
    if (
        type(temperature_K) not in (int, float)
        or not math.isfinite(temperature_K)
        or temperature_K <= 0
    ):
        raise ValueError("temperature_K must be finite and positive")
    if cfg["species"] not in atomic_numbers:
        raise ValueError("Unknown element symbol")
    if timestep_ps is not None and frame_interval_fs is not None:
        raise ValueError("Specify physical timestep or saved-frame interval, not both")
    frames = read(path, ":")
    if len(frames) < 4:
        raise ValueError("At least four frames required")
    tagged = ["physical_step" in a.info for a in frames]
    if any(tagged) and not all(tagged):
        raise ValueError("Partial physical_step metadata")
    if all(tagged):
        if timestep_ps is None or not math.isfinite(timestep_ps) or timestep_ps <= 0:
            raise ValueError("Physical-step trajectories require positive timestep_ps")
        times = (
            np.array([a.info["physical_step"] for a in frames], dtype=float)
            * timestep_ps
        )
        if frame_interval_fs is not None:
            raise ValueError(
                "Specify physical timestep or saved-frame interval, not both"
            )
    else:
        if (
            frame_interval_fs is None
            or not math.isfinite(frame_interval_fs)
            or frame_interval_fs <= 0
        ):
            raise ValueError(
                "Trajectories without physical_step require positive frame_interval_fs"
            )
        times = np.arange(len(frames)) * frame_interval_fs / 1000
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("Frame times must be finite and strictly increasing")
    spacing = np.diff(times)
    if not np.allclose(spacing, spacing[0], rtol=1e-8, atol=1e-12):
        raise ValueError("Uniform saved-frame times required")
    if cfg["smoothed"] == "max" and spacing[0] > 1.0 + 1e-12:
        raise ValueError("pymatgen max smoothing requires saved-frame spacing <= 1 ps")
    selected = times >= cfg["equilibration_ps"] - 1e-10
    frames = [a for a, keep in zip(frames, selected) if keep]
    times = times[selected]
    if len(frames) < 4 or len({a.info.get("system_id") for a in frames}) != 1:
        raise ValueError(
            "One ordered trajectory and at least four production frames required"
        )
    first = frames[0]
    ions = first.numbers == atomic_numbers[cfg["species"]]
    if not ions.any() or ions.all() or first.get_volume() <= 0:
        raise ValueError("Mobile ions and a nonempty framework required")
    if any(
        not a.pbc.all()
        or not np.array_equal(a.numbers, first.numbers)
        or not np.isfinite(a.positions).all()
        or not np.allclose(a.cell, first.cell, rtol=0, atol=1e-9)
        for a in frames
    ):
        raise ValueError(
            "Finite positions, periodic fixed cell and preserved atom order required"
        )
    structures = [AseAtomsAdaptor.get_structure(a) for a in frames]
    fractional = np.array([s.frac_coords for s in structures])
    jumps = np.diff(fractional, axis=0)
    jumps -= np.round(jumps)
    if np.max(np.abs(jumps)) >= 0.45:
        raise ValueError(
            "Saved frame spacing is too large for reliable periodic unwrapping"
        )
    # This guard cannot detect a complete cell traversal between saved frames.
    report = dict(
        schema="atomistic.transport/v1",
        status="provisional",
        species=cfg["species"],
        ionic_charge_e=cfg["charge"],
        temperature_K=temperature_K,
        configuration=cfg,
        equilibration_ps=cfg["equilibration_ps"],
        production_duration_ps=float(times[-1] - times[0]),
        frames=len(frames),
        ion_count=int(ions.sum()),
        volume_angstrom3=float(first.get_volume()),
        provider="pymatgen-analysis-diffusion",
        provider_version=version("pymatgen-analysis-diffusion"),
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        method=(
            "framework-drift-corrected single-origin MSD; unweighted fit"
            if cfg["smoothed"] is False
            else "framework-drift-corrected multi-origin MSD; max smoothing"
        ),
        limitations=[
            "Finite trajectory; diffusion regime and size/time convergence require review.",
            "Fit standard error is not a convergence test; MSD points are correlated.",
            "Nernst-Einstein estimate neglects ion-ion cross correlations.",
        ],
        diffusivity_cm2_s=None,
        conductivity_NE_mS_cm=None,
    )
    try:
        analysis = DiffusionAnalyzer.from_structures(
            structures,
            cfg["species"],
            temperature_K,
            time_step=float(spacing[0] * 1000),
            step_skip=1,
            smoothed=cfg["smoothed"],
            min_obs=cfg["min_observations"],
        )
    except ValueError as error:
        if "Not enough data" not in str(error):
            raise
        return dict(
            report,
            reason="Insufficient production frames for the requested MSD estimator",
        )
    t, msd = analysis.dt, analysis.msd
    mask = np.ones(len(t), dtype=bool)
    if cfg["fit_start_ps"] is not None:
        mask &= t >= cfg["fit_start_ps"] * 1000
    if cfg["fit_end_ps"] is not None:
        mask &= t <= cfg["fit_end_ps"] * 1000
    if mask.sum() < 3:
        raise ValueError("Fit interval must contain at least three lag points")
    ft, fm = t[mask], msd[mask]
    design = np.column_stack([ft, np.ones_like(ft)])
    weights = np.sqrt(1 / ft) if cfg["smoothed"] == "max" else np.ones_like(ft)
    slope, intercept = np.linalg.lstsq(
        design * weights[:, None], fm * weights, rcond=None
    )[0]
    residual = np.sum((fm - (slope * ft + intercept)) ** 2)
    variance = np.sum((fm - fm.mean()) ** 2)
    report.update(
        fit_lag_ps=[float(ft[0] / 1000), float(ft[-1] / 1000)],
        msd_final_angstrom2=float(msd[-1]),
        fit_slope_angstrom2_fs=float(slope),
        fit_r2=float(1 - residual / variance) if variance > 0 else None,
        max_framework_displacement_angstrom=float(analysis.max_framework_displacement),
    )
    if math.isfinite(slope) and slope > 0:
        diffusion, error = get_diffusivity_from_msd(fm, ft, smoothed=cfg["smoothed"])
        factor = (
            int(ions.sum())
            / (first.get_volume() * 1e-30)
            * (cfg["charge"] * elementary_charge) ** 2
            / (Boltzmann * temperature_K)
        )
        report.update(
            diffusivity_cm2_s=float(diffusion),
            conductivity_NE_mS_cm=float(diffusion * 1e-4 * factor * 10),
            fit_standard_error_cm2_s=float(error),
        )
    else:
        report["reason"] = (
            "Nonpositive fitted slope; no established diffusion coefficient"
        )
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        np.savetxt(
            out / "msd.csv",
            np.column_stack([t / 1000, msd]),
            delimiter=",",
            header="lag_ps,msd_angstrom2",
            comments="",
        )
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        with plt.rc_context({"font.size": 14}):
            fig, axis = plt.subplots(figsize=(6, 5))
            axis.plot(t / 1000, msd, label=cfg["species"] + " MSD", linewidth=2.5)
            axis.plot(
                ft / 1000, slope * ft + intercept, "--", label="fit", linewidth=2.5
            )
            axis.set_xlabel("Lag time (ps)", fontweight="bold")
            axis.set_ylabel("MSD (angstrom²)", fontweight="bold")
            axis.set_title(f"{temperature_K:g} K — provisional")
            axis.legend(frameon=False)
            axis.grid(False)
            fig.tight_layout()
            for suffix in ("png", "svg"):
                fig.savefig(out / f"msd.{suffix}", dpi=160, bbox_inches="tight")
            plt.close(fig)
    return report
