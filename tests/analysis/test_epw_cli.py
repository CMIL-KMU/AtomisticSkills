"""Shared table parsing preserves each standalone scientific CLI's output."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / ".agents/skills/mat-epw-mobility/scripts"


@pytest.mark.parametrize("script", ["parse_prt.py", "parse_epw_prtgkk.py"])
def test_standalone_table_parser_from_other_directory(tmp_path, script):
    source = tmp_path / "phonon.out"
    source.write_text(
        "q coord.: 0.1 0.2 0.3\niq = 1 coord.: 0.1 0.2 0.3\n"
        "k coord.: 0.0 0.0 0.0\n"
        "2 2 1 -3 -3 12 4\n2 2 2 -3 -3 15 5\n"
        "3 3 1 10 10 12 99\nnoise\n"
    )
    process = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / script), str(source), "--fermi", "-3"],
        cwd=tmp_path, check=True, text=True, capture_output=True,
    )
    result = json.loads(process.stdout)
    assert result["q_cartesian"] == [0.1, 0.2, 0.3]
    assert result["n_rows"] == 3
    assert result["cbm_ibnd"] == 2
    assert result["rank1_g_meV"] == 5
    assert result["sum_g2_LO_TO_quartet_meV2"] == 41
