"""Plan or install source-declared dependencies into new, explicit Conda prefixes.

This does not remove named environments or execute the legacy install.sh scripts.
Package installation is distinct from engine/checkpoint and scientific qualification.
Run in an environment containing PyYAML. Logs and receipts stay under --directory.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import yaml

SOURCE = Path(__file__).resolve().parents[1]


def plan(name):
    """Read the actual environment declaration without guessing missing providers."""
    root = SOURCE / "conda-envs" / name
    declaration = root / "core_env.yaml"
    if not declaration.exists():
        declaration = root / "env.yaml"
    data = yaml.safe_load(declaration.read_text())
    pip = [item for group in data["dependencies"] if isinstance(group, dict)
           for item in group.get("pip", [])]
    missing = [item for item in pip if item.startswith("-e /")
               and not Path(item[3:].split("[")[0]).exists()]
    return dict(name=name, declaration=str(declaration.relative_to(SOURCE)),
                source_sha256=hashlib.sha256(declaration.read_bytes()).hexdigest(),
                conda={k: v for k, v in data.items() if k not in ("name", "prefix")},
                pip=pip, missing_local_sources=missing,
                instructions=str((root / "README.md").relative_to(SOURCE)),
                qualification="not-tested")


def install(item, directory, conda):
    """Install dependencies once; never reuse an unreceipted or different prefix."""
    work = directory / item["name"]
    try:
        work.mkdir()
    except FileExistsError:
        return dict(name=item["name"], status="existing-prefix-or-receipt",
                    detail="Use a new directory; no existing installation was changed")
    receipt = work / "receipt.json"
    prefix = work / "environment"
    state = dict(name=item["name"], source_sha256=item["source_sha256"],
                 prefix=str(prefix), status="pending", qualification="not-tested")

    def save():
        state["updated_at"] = time.time()
        temporary = receipt.with_suffix(".tmp")
        temporary.write_text(json.dumps(state, indent=2) + "\n")
        temporary.replace(receipt)

    if item["missing_local_sources"]:
        state.update(status="blocked-local-source", missing=item["missing_local_sources"])
        save()
        return state
    configuration = dict(item["conda"])
    configuration["dependencies"] = [d for d in configuration["dependencies"] if isinstance(d, str)]
    declaration = work / "conda.yaml"
    declaration.write_text(yaml.safe_dump(configuration, sort_keys=False))
    requirements = work / "requirements.txt"
    requirements.write_text("\n".join(item["pip"]) + "\n")
    env = dict(os.environ, CONDA_PKGS_DIRS=str(directory / "package-cache"),
               PIP_CACHE_DIR=str(directory / "pip-cache"), PYTHONNOUSERSITE="1",
               PATH=str(prefix / "bin") + os.pathsep + os.environ.get("PATH", ""))
    commands = [[conda, "env", "create", "--prefix", str(prefix), "--file", str(declaration), "--yes"]]
    if item["pip"]:
        commands.append([str(prefix / "bin/python"), "-m", "pip", "install", "-r", str(requirements)])
    for ordinal, command in enumerate(commands):
        state.update(status="installing", stage=ordinal)
        save()
        with (work / f"install-{ordinal}.log").open("w") as log:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            state.update(status="installation-failed", exit_code=result.returncode)
            save()
            return state
    check = subprocess.run([str(prefix / "bin/python"), "-m", "pip", "check"],
                           env=env, capture_output=True, text=True)
    (work / "pip-check.log").write_text(check.stdout + check.stderr)
    freeze = subprocess.run([str(prefix / "bin/python"), "-m", "pip", "freeze"],
                            env=env, capture_output=True, text=True, check=True)
    (work / "freeze.txt").write_text(freeze.stdout)
    state.update(status="dependencies-installed" if check.returncode == 0 else "dependency-conflicts",
                 qualification="engine-model-and-skill-tests-required",
                 postinstall_instructions=item["instructions"])
    save()
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "install"))
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--conda", default="conda")
    parser.add_argument("--workers", type=int, choices=(1, 2), default=2)
    parser.add_argument("--environments", nargs="+")
    args = parser.parse_args()
    directory = args.directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    available = {p.name for p in (SOURCE / "conda-envs").iterdir() if p.is_dir()}
    names = args.environments or sorted(available)
    if not set(names) <= available or len(names) != len(set(names)):
        parser.error("Unknown or duplicate source environment")
    planned = [plan(name) for name in names]
    (directory / "plan.json").write_text(json.dumps(planned, indent=2) + "\n")
    if args.action == "install":
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(lambda item: install(item, directory, args.conda), planned))
        (directory / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
        return 0 if all(r["status"] == "dependencies-installed" for r in results) else 1
    else:
        print(json.dumps(planned, indent=2))


if __name__ == "__main__":
    sys.exit(main())
