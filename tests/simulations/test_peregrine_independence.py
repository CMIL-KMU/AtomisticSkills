"""Use a real provider registry/reference fixture while forbidding consumer imports."""
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest


def test_real_cli_without_workflow_imports(tmp_path):
    location = os.environ.get("ATOMISTIC_PEREGRINE_FIXTURE")
    if not location:
        pytest.skip("Explicit real Peregrine registry and reference fixture required")
    fixture_root = Path(location).resolve()
    fixture = json.loads((fixture_root / "fixture.json").read_text())
    source = Path(__file__).resolve().parents[2]
    output = tmp_path / "simulation"
    request = dict(operation="static", structures=fixture["structures"], options={},
                   seed=91, model=fixture["model"])
    payload = dict(request=request, registry=str(fixture_root / "registry"), output=str(output))
    guard = """import importlib.abc,runpy,sys
class IndependentOnly(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0].startswith(('anvil','mkite')):
            raise AssertionError('Consumer import forbidden: '+fullname)
sys.meta_path.insert(0,IndependentOnly())
sys.argv=sys.argv[1:]
runpy.run_path(sys.argv[0],run_name='__main__')
"""
    result = subprocess.run([sys.executable, "-I", "-B", "-c", guard,
        str(source / "tools/run_simulation.py"), "peregrine-run"],
        input=json.dumps(payload), text=True, capture_output=True, timeout=180,
        env=dict(os.environ, NUMBA_CACHE_DIR=str(tmp_path / "numba-cache")))
    assert result.returncode == 0, result.stderr
    observed = json.loads(result.stdout)
    assert observed["schema"] == "atomistic-skills.simulation/v1"
    assert not {"job", "attempt", "runtime", "model_version"} & observed["request"].keys()
    for actual, expected in zip(observed["results"], fixture["expected"], strict=True):
        for key in ("energy", "forces", "stress"):
            np.testing.assert_allclose(actual[key], expected[key], rtol=1e-10, atol=1e-10)
    assert json.loads((output / "simulation.json").read_text()) == observed
    assert not (output / "manifest.json").exists()
