import math
from copy import deepcopy
import pytest
from atomistic_analysis.arrhenius import fit_transport, KB_EV_K


def observations():
    return [
        dict(
            source_uuid=str(t),
            temperature_K=t,
            diffusivity_cm2_s=0.01 * math.exp(-0.3 / (KB_EV_K * t)),
            conductivity_NE_mS_cm=4e5 / t * math.exp(-0.3 / (KB_EV_K * t)),
        )
        for t in (600, 800, 1000, 1200)
    ]


def test_known_barrier_and_conductivity_temperature_factor():
    result = fit_transport(observations())
    assert result["available"]
    for fit in result["fits"]:
        assert fit["activation_energy_eV"] == pytest.approx(0.3)
        assert fit["standard_error_eV"] < 1e-14 and fit["r2"] == pytest.approx(1)
    assert result["fits"][0]["log10_prefactor"] == pytest.approx(-2)
    assert result["fits"][1]["log10_prefactor"] == pytest.approx(math.log10(4e5))


@pytest.mark.parametrize("value", [0, -1, None, float("nan"), float("inf"), True])
def test_invalid_observations_are_excluded_with_reasons(value):
    rows = observations()
    rows[0]["diffusivity_cm2_s"] = value
    result = fit_transport(rows)
    assert result["available"] and len(result["excluded"]) == 1
    assert len(result["accepted_sources"]) == 3


def test_no_duplicate_trajectory_temperature_or_two_point_certainty():
    rows = observations()
    assert not fit_transport(rows[:2])["available"]
    other = deepcopy(rows)
    other[0]["source_uuid"] = other[1]["source_uuid"]
    assert "twice" in fit_transport(other)["reason"]
    rows[0]["temperature_K"] = rows[1]["temperature_K"]
    assert "Repeated" in fit_transport(rows)["reason"]


def test_disordered_data_fit_and_negative_barrier_warning():
    rows = observations()
    assert fit_transport(rows[::-1])["fits"][0][
        "activation_energy_eV"
    ] == pytest.approx(0.3)
    for r in rows:
        r["diffusivity_cm2_s"] = 1 / r["diffusivity_cm2_s"]
    assert fit_transport(rows)["fits"][0]["warning"]


def test_optional_target_and_optional_conductivity():
    rows = observations()
    result = fit_transport(rows)
    assert all("prediction" not in f for f in result["fits"])
    result = fit_transport(rows, target_temperature_K=298.15)
    assert result["fits"][0]["prediction"]["extrapolated"]
    expected = math.log10(0.01) - 0.3 / (KB_EV_K * 298.15 * math.log(10))
    assert result["fits"][0]["prediction"]["log10_value"] == pytest.approx(expected)
    for r in rows:
        r.pop("conductivity_NE_mS_cm")
    result = fit_transport(rows)
    assert result["available"] and len(result["fits"]) == 1
    assert result["unavailable_quantities"][0]["excluded"]
    with pytest.raises(ValueError):
        fit_transport(rows, target_temperature_K=0)


def test_no_fit_after_exclusions_and_no_mixed_species():
    rows = observations()
    rows[0]["diffusivity_cm2_s"] = 0
    rows[1]["diffusivity_cm2_s"] = -1
    assert not fit_transport(rows)["available"]
    rows = observations()
    for r in rows:
        r["species"] = "Na"
    rows[0]["species"] = "Li"
    assert "species" in fit_transport(rows)["reason"]
