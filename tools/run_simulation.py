"""Independent VASP/Peregrine JSON CLI for humans, agents and workflow consumers.

Usage: python tools/run_simulation.py ACTION --input request.json --output result.json
Input/output default to stdin/stdout. Run in the selected provider environment.
"""
import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def dispatch(action: str, data: dict):
    """Dispatch explicit scientific operations without consumer credentials or IDs."""
    if action == "identity":
        from src.simulations.identity import identity
        return identity(data["provider"])
    if action == "peregrine-run":
        from src.simulations.peregrine.runner import execute
        return execute(data["request"], data["registry"], data["output"])
    if action == "peregrine-model-precision":
        from src.simulations.peregrine.runner import same_at_model_precision
        return bool(same_at_model_precision(**data))
    if action.startswith("vasp-"):
        from src.simulations.vasp import protocols, dynamics
        data = dict(data)
        directory = data.pop("protocol_directory", None)
        if directory is not None:
            protocols.ROOT = Path(directory).resolve()
        if action == "vasp-registry":
            return protocols.registry()
        if action == "vasp-prepare":
            return protocols.prepare(**data)
        if action == "vasp-parse":
            if "stage" in data:
                return dynamics.parse(**data)
            from src.simulations.vasp.parser import parse_stage
            return parse_stage(**data)
        if action == "vasp-setup":
            dynamics.setup(Path(data["folder"]), data["stage"])
            return True
        if action == "vasp-inputs":
            from src.simulations.vasp.inputs import write_inputs
            return write_inputs(**data)
        if action == "vasp-neb-observables":
            return dynamics.neb_observables(data["path"])
        if action == "vasp-mobile-projection":
            from pymatgen.core import Structure
            return dynamics.mobile_projection(data["values"], Structure.from_dict(data["structure"])).tolist()
        if action == "vasp-ldau":
            from src.simulations.vasp.ldau import apply
            return apply(data["prepared"], data["selection"])
    raise ValueError("Unknown simulation action")


def main() -> None:
    """Keep Python/native library diagnostics off the JSON response channel."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("identity", "vasp-registry", "vasp-prepare", "vasp-parse",
        "vasp-setup", "vasp-inputs", "vasp-neb-observables", "vasp-mobile-projection", "vasp-ldau",
        "peregrine-run", "peregrine-model-precision"))
    parser.add_argument("--input", type=Path, help="JSON request; defaults to stdin")
    parser.add_argument("--output", type=Path, help="JSON response; defaults to stdout")
    args = parser.parse_args()
    data = json.loads(args.input.read_text()) if args.input else json.load(sys.stdin)
    with os.fdopen(os.dup(sys.stdout.fileno()), "w") as response:
        sys.stdout.flush()
        os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
        sys.stdout = sys.stderr
        result = dispatch(args.action, data)
        from monty.json import MontyEncoder
        text = json.dumps(result, cls=MontyEncoder, allow_nan=False)
        if args.output:
            with args.output.open("x") as output:
                output.write(text + "\n")
        else:
            response.write(text + "\n")


if __name__ == "__main__":
    main()
