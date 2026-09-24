"""Measure independent simulation source and numerical dependency versions."""
import hashlib
import json
from importlib.metadata import version
from pathlib import Path


def identity(provider: str) -> dict:
    """Return bytes and dependencies that determine this provider's behavior."""
    if provider not in {"peregrine"}:
        raise ValueError("Unknown simulation provider")
    root = Path(__file__).resolve().parents[2]
    paths = [Path(__file__), root / "tools/run_simulation.py"]
    paths += sorted((Path(__file__).parent / provider).glob("*.py"))
    paths += sorted((Path(__file__).parent / provider).glob("*.json"))
    if provider == "peregrine":
        paths += sorted((root / "src/utils/analysis").glob("*.py"))
    files = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    dependencies = (
        "ase", "numpy", "torch", "peregrine-pot", "peregrine-sim", "jsonschema")
    return dict(files=files, sha256=hashlib.sha256(json.dumps(files,sort_keys=True).encode()).hexdigest(),
                versions={name:version(name) for name in dependencies})
