"""Original script entry points and the neutral CLI run without a project wheel."""

import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_script_from_an_unrelated_working_directory(tmp_path):
    script = ROOT / ".agents/skills/mat-diffusion-analysis/scripts/calculate_activation_energy.py"
    rows = [dict(source_uuid=str(t), temperature_K=t, diffusivity_cm2_s=d)
            for t, d in [(800, 1e-5), (1000, 2e-5), (1200, 3e-5)]]
    source = tmp_path / "observations.json"
    source.write_text(json.dumps(rows))
    output = tmp_path / "result"
    subprocess.run([sys.executable, "-I", "-B", str(script), "--observations", str(source),
                    "--output_dir", str(output)], cwd=tmp_path, check=True, capture_output=True)
    report = json.loads((output / "arrhenius.json").read_text())
    assert report["available"] and len(report["fits"]) == 1
    assert (output / "arrhenius-0.png").is_file()
    assert json.loads(source.read_text()) == rows


def test_neutral_cli_keeps_drug_tools_and_rejects_arbitrary_dispatch(tmp_path):
    command = [sys.executable, "-I", "-B", str(ROOT / "tools/run_analysis.py")]
    result = subprocess.run([*command, "describe"], input="{}", text=True,
                            cwd=tmp_path, capture_output=True, check=True)
    descriptor = json.loads(result.stdout)
    assert "docking-box" in descriptor["schemas"]
    assert "supercell" in descriptor["schemas"]
    invalid = subprocess.run([*command, "import-os"], input="{}", text=True,
                             cwd=tmp_path, capture_output=True)
    assert invalid.returncode != 0


def test_source_identity_covers_runner_and_shared_science():
    from src.utils.analysis.tools import source_identity
    identity = source_identity()
    assert "tools/run_analysis.py" in identity["files"]
    assert "src/utils/analysis/transport.py" in identity["files"]
    assert "src/utils/analysis/tools.py" in identity["files"]
    assert len(identity["sha256"]) == 64
