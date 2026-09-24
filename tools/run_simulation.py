"""Independent Peregrine JSON CLI for humans, agents and workflow consumers.

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
    raise ValueError("Unknown simulation action")


def main() -> None:
    """Keep Python/native library diagnostics off the JSON response channel."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("identity", "peregrine-run", "peregrine-model-precision"))
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
