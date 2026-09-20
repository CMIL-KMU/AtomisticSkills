"""Arrhenius diagnostics of explicit observations; no database or folder discovery."""

import hashlib
import math
from pathlib import Path
from statistics import linear_regression, mean

KB_EV_K = 1.380649e-23 / 1.602176634e-19


def positive(value):
    """Exclude booleans, missing data and nonfinite/nonpositive numbers."""
    return type(value) in (int, float) and math.isfinite(value) and value > 0


def fit_transport(rows: list[dict], *, target_temperature_K=None) -> dict:
    """Fit log D and log(sigma_NE T); disclose exclusions and optional extrapolation.

    Each row needs source_uuid, temperature_K and diffusivity_cm2_s. Conductivity is
    optional. The caller selects comparable material/method observations explicitly.
    At least three independent temperatures are needed for regression scatter.
    """
    if target_temperature_K is not None and not positive(target_temperature_K):
        raise ValueError("target_temperature_K must be finite and positive")
    result = dict(
        schema="atomistic.arrhenius/v1",
        status="provisional",
        available=False,
        fits=[],
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        method="Unweighted least squares in log10 space; one observation per temperature",
        selected_count=len(rows),
        excluded=[],
        target_temperature_K=target_temperature_K,
    )
    ids = [r.get("source_uuid") for r in rows]
    if any(not isinstance(i, str) or not i for i in ids):
        return dict(result, reason="Every observation needs a source_uuid")
    if len(set(ids)) != len(ids):
        return dict(result, reason="The same MD trajectory cannot be counted twice.")
    for key in ("species", "ionic_charge_e", "method"):
        values = [r.get(key) for r in rows]
        if any(v is not None for v in values) and len(set(values)) > 1:
            return dict(
                result, reason=f"Incompatible or incomplete {key} across observations"
            )
    accepted = []
    for row in rows:
        bad = next(
            (
                k
                for k in ("temperature_K", "diffusivity_cm2_s")
                if not positive(row.get(k))
            ),
            None,
        )
        if bad:
            result["excluded"].append(
                dict(
                    source_uuid=row["source_uuid"],
                    temperature_K=(
                        row.get("temperature_K")
                        if positive(row.get("temperature_K"))
                        else None
                    ),
                    field=bad,
                    reason="Missing, nonfinite or nonpositive value; excluded from log fit",
                )
            )
        else:
            accepted.append(row)
    result["accepted_sources"] = [r["source_uuid"] for r in accepted]
    if len(accepted) < 3:
        return dict(
            result,
            reason="At least three eligible distinct temperatures required after exclusions.",
        )
    temps = [r["temperature_K"] for r in accepted]
    if len(set(temps)) != len(temps):
        return dict(
            result,
            reason="Repeated temperatures require explicit replicate handling; select one analysis per temperature.",
        )
    result["temperature_range_K"] = [min(temps), max(temps)]
    for name, field, unit, formula, scale_t in (
        ("Diffusivity", "diffusivity_cm2_s", "cm²/s", "D = D0 exp(-Ea / kBT)", False),
        (
            "Nernst–Einstein conductivity × T",
            "conductivity_NE_mS_cm",
            "mS K/cm",
            "σNE T = A exp(-Ea / kBT)",
            True,
        ),
    ):
        points, excluded = [], []
        for row in accepted:
            value = row.get(field)
            value = (
                value * row["temperature_K"] if scale_t and positive(value) else value
            )
            if not positive(value):
                excluded.append(
                    dict(
                        source_uuid=row["source_uuid"],
                        field=field,
                        reason="Missing, nonfinite or nonpositive value",
                    )
                )
            else:
                points.append((row, value))
        if len(points) < 3:
            result.setdefault("unavailable_quantities", []).append(
                dict(
                    name=name,
                    reason="Fewer than three eligible temperatures",
                    excluded=excluded,
                )
            )
            continue
        x = [1000 / r["temperature_K"] for r, _ in points]
        y = [math.log10(v) for _, v in points]
        slope, intercept = linear_regression(x, y)
        predicted = [slope * v + intercept for v in x]
        sse = sum((a - b) ** 2 for a, b in zip(y, predicted))
        sst = sum((v - mean(y)) ** 2 for v in y)
        sxx = sum((v - mean(x)) ** 2 for v in x)
        factor = 1000 * KB_EV_K * math.log(10)
        fit = dict(
            name=name,
            formula=formula,
            unit=unit,
            activation_energy_eV=-slope * factor,
            standard_error_eV=math.sqrt(sse / (len(x) - 2) / sxx) * factor,
            r2=1 - sse / sst if sst else None,
            log10_prefactor=intercept,
            slope=slope,
            intercept=intercept,
            points=sorted(zip(x, y)),
            line=[[v, slope * v + intercept] for v in sorted((min(x), max(x)))],
            accepted_sources=[r["source_uuid"] for r, _ in points],
            excluded=excluded,
            warning="Nonpositive activation energy; inspect selected observations."
            if slope >= 0
            else None,
        )
        if target_temperature_K is not None:
            value = slope * 1000 / target_temperature_K + intercept
            if scale_t:
                value -= math.log10(target_temperature_K)
            fit["prediction"] = dict(
                temperature_K=target_temperature_K,
                log10_value=value,
                unit="mS/cm" if scale_t else "cm²/s",
                extrapolated=not min(r["temperature_K"] for r, _ in points)
                <= target_temperature_K
                <= max(r["temperature_K"] for r, _ in points),
                warning="Model prediction, not a measurement; assumes the same transport mechanism.",
            )
        result["fits"].append(fit)
    return dict(result, available=bool(result["fits"]))
