"""CLI adapters for explicit scientific inputs; never submit jobs or mutate a DB."""

import argparse
import hashlib
import json
import shutil
from tempfile import TemporaryDirectory
from pathlib import Path


def save(directory, name, result, settings):
    """Persist effective inputs separately from output, without replacing earlier runs."""
    import yaml

    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    (out / "input_configs.yaml").write_text(yaml.safe_dump(settings, sort_keys=True))
    (out / name).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


def diffusion_main():
    """Analyze one trajectory with user-selected ion, charge and explicit timing."""
    from src.utils.analysis.transport import analyze, configuration

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trajectory")
    parser.add_argument("--species", required=True)
    parser.add_argument("--charge", type=int, required=True)
    parser.add_argument("--temperature", type=float, required=True)
    timing = parser.add_mutually_exclusive_group(required=True)
    timing.add_argument("--frame-interval-fs", type=float)
    timing.add_argument("--timestep-ps", type=float)
    parser.add_argument("--ignore_ps", type=float, required=True)
    parser.add_argument("--smoothed", choices=["false", "max"], default="false")
    parser.add_argument("--min-observations", type=int, default=30)
    parser.add_argument("--fit-start-ps", type=float)
    parser.add_argument("--fit-end-ps", type=float)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()
    if Path(args.output_dir).exists():
        parser.error("Use a new output directory; prior analyses are immutable")
    cfg = configuration(
        dict(
            species=args.species,
            charge=args.charge,
            equilibration_ps=args.ignore_ps,
            smoothed=False if args.smoothed == "false" else "max",
            min_observations=args.min_observations,
            fit_start_ps=args.fit_start_ps,
            fit_end_ps=args.fit_end_ps,
        )
    )
    with TemporaryDirectory() as temporary:
        result = analyze(
            args.trajectory,
            cfg,
            temperature_K=args.temperature,
            timestep_ps=args.timestep_ps,
            frame_interval_fs=args.frame_interval_fs,
            output_dir=temporary,
        )
        result["source_sha256"] = hashlib.sha256(
            Path(args.trajectory).read_bytes()
        ).hexdigest()
        save(
            args.output_dir,
            "diffusion_results.json",
            result,
            dict(vars(args), effective_configuration=cfg),
        )
        for item in Path(temporary).iterdir():
            shutil.copyfile(item, Path(args.output_dir) / item.name)


def arrhenius_main():
    """Fit explicitly listed observations; extrapolation is opt-in."""
    from src.utils.analysis.arrhenius import fit_transport

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--observations",
        required=True,
        help="JSON list with source_uuid, temperature_K and diffusivity_cm2_s; optional conductivity_NE_mS_cm",
    )
    parser.add_argument("--target-temperature-K", type=float)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()
    rows = json.loads(Path(args.observations).read_text())
    result = fit_transport(rows, target_temperature_K=args.target_temperature_K)
    result["source_sha256"] = hashlib.sha256(
        Path(args.observations).read_bytes()
    ).hexdigest()
    save(args.output_dir, "arrhenius.json", result, vars(args))
    if not result["available"]:
        parser.exit(2, result["reason"] + "\n")
    plot_arrhenius(result, args.output_dir)


def plot_arrhenius(result, directory):
    """Render the existing fit points/lines without recomputing scientific values."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from src.utils.analysis.plotting import plot_style

    with plot_style():
        for index, fit in enumerate(result["fits"]):
            fig, axis = plt.subplots(figsize=(6, 5))
            axis.scatter(*zip(*fit["points"]), label="accepted observations")
            axis.plot(
                *zip(*fit["line"]),
                linewidth=2.5,
                marker="",
                color="#ee9b00",
                label=f"Ea = {fit['activation_energy_eV']:.3f} eV",
            )
            axis.set_xlabel("1000 / T (K⁻¹)", fontweight="bold")
            axis.set_ylabel(f"log10 [{fit['unit']}]", fontweight="bold")
            axis.set_title(f"{fit['name']} — provisional", fontsize=11)
            axis.legend(frameon=False, fontsize=11)
            axis.grid(False)
            fig.tight_layout()
            for suffix in ("png", "svg"):
                fig.savefig(
                    Path(directory) / f"arrhenius-{index}.{suffix}",
                    dpi=160,
                    bbox_inches="tight",
                )
            plt.close(fig)
