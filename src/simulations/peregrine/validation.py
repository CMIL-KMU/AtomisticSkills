"""Portable request validation shared by database and scientific workers."""
import json
from pathlib import Path
from jsonschema import Draft202012Validator, FormatChecker

WORKER_SCHEMA = json.loads((Path(__file__).with_name("request.json")).read_text())


def validate_request_api(request):
    """Reject removed API choices before persistence, preflight or execution."""
    if "engine" in request or "engine" in request.get("options", {}):
        raise ValueError("engine was removed; omit it to use torch-sim. ASE simulation/restart is no longer supported")
    if "properties" in request or "properties" in request.get("options", {}):
        raise ValueError("properties is not a simulation option; The runner configures energy/forces/stress through PeregrineModel")
    if "restart" in request and set(request["restart"]) != {"manifest_sha256"}:
        raise ValueError("Unsupported restart request")
    json.dumps(request, allow_nan=False)
    if "operation" in request:
        errors = list(Draft202012Validator(WORKER_SCHEMA, format_checker=FormatChecker()).iter_errors(request))
        if errors:
            raise ValueError("Invalid worker operation: " + errors[0].message)
    if "transport" in request and (request.get("operation") != "md" or len(request["structures"]) != 1 or "bias" in request):
        raise ValueError("Transport requires one unbiased MD trajectory")
    if "checkpoint_interval" in request and request.get("operation") != "md":
        raise ValueError("Periodic checkpoints require MD")
    if "bias" in request:
        validate_bias(request)


def validate_bias(request):
    """Check dimensions that cannot be expressed by the shared JSON Schema."""
    bias = request['bias']
    errors = list(Draft202012Validator(WORKER_SCHEMA['properties']['bias']).iter_errors(bias))
    if errors:
        raise ValueError('Invalid bias: ' + errors[0].message)
    structures, cvs = request['structures'], bias['cvs']
    if request.get('operation') != 'md' or len(request.get('system_ids', [])) != len(structures):
        raise ValueError('Bias requires MD and ordered system IDs')
    for cv in cvs:
        if any(max(cv['atoms']) >= len(row['numbers']) for row in structures):
            raise ValueError('CV local indices must exist in every structure')
        if cv['kind'] == 'position_projection' and not any(cv['direction']):
            raise ValueError('Projection direction must be nonzero')
    matrices = [bias[k] for k in ('centers', 'force_constants')] if bias['kind'] == 'umbrella' else [[bias['sigma']]]
    for matrix in matrices:
        if len(matrix) != (len(structures) if bias['kind'] == 'umbrella' else 1) or any(len(row) != len(cvs) for row in matrix):
            raise ValueError('Bias arrays require explicit system-by-CV dimensions')


