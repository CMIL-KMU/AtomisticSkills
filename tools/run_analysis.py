"""Call existing scientific helpers with explicit JSON inputs, without installation.

This script is a neutral CLI for external consumers and standalone users. It
does not know any database, campaign, scheduler, credentials or consumer policy.
"""

import argparse
from contextlib import redirect_stdout
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.analysis import tools


def dispatch(action, data):
    if action == "describe":
        return dict(schemas=tools.SCHEMAS, dependencies=tools.DEPENDENCIES)
    if action == "source":
        return tools.source_identity()
    if action == "identity":
        return tools.identity(data["tool"])
    if action == "dependency":
        return tools.dependency_version(data["name"])
    if action == "dependencies":
        from importlib.metadata import PackageNotFoundError
        versions = {}
        names = {"jsonschema", "ase", "numpy", "scipy", "pymatgen", "pymatgen-analysis-diffusion"}
        names.update(name for group in tools.DEPENDENCIES.values() for name in group)
        for name in sorted(names):
            try:
                versions[name] = tools.dependency_version(name)
            except PackageNotFoundError:
                versions[name] = None
        return versions
    if action == "validate":
        return tools.validate(data["tool"], data["inputs"])
    if action == "validate-report":
        tools.validate_report(data["tool"], data["report"])
        return True
    if action == "run":
        return tools.run(data["tool"], data["inputs"])
    if action == "transport-config":
        from src.utils.analysis.transport import configuration
        return configuration(data)
    if action == "transport":
        from src.utils.analysis.transport import analyze
        return analyze(**data)
    if action == "arrhenius":
        from src.utils.analysis.arrhenius import fit_transport
        return fit_transport(**data)
    if action == "plot-arrhenius":
        from src.utils.analysis.cli import plot_arrhenius
        plot_arrhenius(data["result"], data["directory"])
        return True
    if action == "browser-style":
        from src.utils.analysis.plotting import browser_style
        return browser_style()
    raise ValueError("Unknown scientific action")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("describe", "source", "identity", "dependency", "dependencies",
        "validate", "validate-report", "run", "transport-config", "transport",
        "arrhenius", "plot-arrhenius", "browser-style", "structure-png"))
    args = parser.parse_args()
    if args.action == "structure-png":
        from src.utils.analysis.structure_cli import main as render
        render()
        return
    data = json.load(sys.stdin)
    # Library progress output must not corrupt the machine-readable result.
    with redirect_stdout(sys.stderr):
        result = dispatch(args.action, data)
    json.dump(result, sys.stdout, allow_nan=False)


if __name__ == "__main__":
    main()
