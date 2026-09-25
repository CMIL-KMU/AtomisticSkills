"""CLI outputs preserve source files, effective configuration and opt-in prediction."""

import json
import math
import sys
import pytest


def test_diffusion_cli_and_no_overwrite(tmp_path, monkeypatch):
    from ase import Atoms
    from ase.io import write
    from src.utils.analysis.cli import diffusion_main

    path = tmp_path / "source.xyz"
    frames = [
        Atoms(
            "NaCl",
            positions=[[math.sqrt(i) * 0.02, 0, 0], [5, 5, 5]],
            cell=[20] * 3,
            pbc=True,
        )
        for i in range(100)
    ]
    write(path, frames)
    original = path.read_bytes()
    output = tmp_path / "result"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "diffusion",
            str(path),
            "--species",
            "Na",
            "--charge",
            "1",
            "--temperature",
            "700",
            "--frame-interval-fs",
            "10",
            "--ignore_ps",
            "0",
            "--output_dir",
            str(output),
        ],
    )
    diffusion_main()
    report = json.loads((output / "diffusion_results.json").read_text())
    assert report["configuration"]["smoothed"] is False and report["species"] == "Na"
    assert (output / "msd.png").stat().st_size > 0 and (
        output / "msd.svg"
    ).stat().st_size > 0
    assert (output / "input_configs.yaml").is_file()
    assert path.read_bytes() == original
    with pytest.raises(SystemExit):
        diffusion_main()


def test_arrhenius_cli_uses_explicit_observations(tmp_path, monkeypatch):
    from src.utils.analysis.cli import arrhenius_main

    path = tmp_path / "observations.json"
    out = tmp_path / "fit"
    path.write_text(
        json.dumps(
            [
                dict(source_uuid=str(t), temperature_K=t, diffusivity_cm2_s=d)
                for t, d in [(600, 0), (800, 1e-5), (1000, 2e-5), (1200, 3e-5)]
            ]
        )
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "arrhenius",
            "--observations",
            str(path),
            "--target-temperature-K",
            "298.15",
            "--output_dir",
            str(out),
        ],
    )
    arrhenius_main()
    report = json.loads((out / "arrhenius.json").read_text())
    assert report["available"] and len(report["excluded"]) == 1
    assert report["fits"][0]["prediction"]["temperature_K"] == 298.15
    assert (out / "arrhenius-0.png").is_file() and (
        out / "input_configs.yaml"
    ).is_file()
