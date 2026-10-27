"""
Property tests for the liquefaction models. The probability that the factor of safety is below
fsl should increase with fsl, peak acceleration, and magnitude, and decrease with penetration
resistance. A known-value test for Boulanger and Idriss (2016) is in test_hazard_integration.py.
"""

import copy

import numpy as np
import pytest

from ucla_plha import plha

MODELS = {
    "cetin_et_al_2018": (
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "vs12": 180.0,
            "depth": 5.0,
            "n160": 12.0,
            "fc": 8.0,
            "pa": 101.325,
        },
        "n160",
    ),
    "moss_et_al_2006": (
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "depth": 5.0,
            "qc": 5.0,
            "fs": 0.05,
            "pa": 101.325,
        },
        "qc",
    ),
    "boulanger_idriss_2012": (
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "depth": 5.0,
            "n160": 14.0,
            "fc": 5.0,
            "pa": 101.325,
        },
        "n160",
    ),
    "boulanger_idriss_2016": (
        {"sigmav": 100.0, "sigmavp": 60.0, "depth": 5.0, "qc1ncs": 90.0, "pa": 101.325},
        "qc1ncs",
    ),
    "ngl_smt_2024": (
        {
            "ztop": [1.0, 3.0, 6.0],
            "zbot": [3.0, 6.0, 9.0],
            "qc1ncs": [60.0, 80.0, 100.0],
            "ic": [1.8, 1.9, 2.0],
            "sigmav": [38.0, 85.0, 140.0],
            "sigmavp": [28.0, 55.0, 85.0],
            "ksat": [1.0, 1.0, 1.0],
            "pa": 101.325,
        },
        "qc1ncs",
    ),
}
FSL = np.array([0.5, 0.8, 1.0, 1.25, 1.6])
PGA = np.array([0.1, 0.2, 0.3, 0.45])


def _cdfs(model, params, pga, m=7.0, sigma=0.5):
    n = len(pga)
    config = {"liquefaction_models": {model: params}}
    cdfs, _ = plha.get_liquefaction_cdfs(
        np.full(n, m), np.log(pga), np.full(n, sigma), FSL, model, config
    )
    return cdfs


def _scale(value, factor):
    if isinstance(value, list):
        return [v * factor for v in value]
    return value * factor


@pytest.mark.parametrize("model", sorted(MODELS))
def test_cdf_increases_with_fsl(model):
    params, _ = MODELS[model]
    cdfs = _cdfs(model, params, PGA)
    assert cdfs.shape == (len(PGA), len(FSL))
    assert np.all(np.diff(cdfs, axis=1) > 0.0)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_cdf_increases_with_pga(model):
    params, _ = MODELS[model]
    cdfs = _cdfs(model, params, PGA)
    assert np.all(np.diff(cdfs, axis=0) > 0.0)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_cdf_increases_with_magnitude(model):
    params, _ = MODELS[model]
    small = _cdfs(model, params, PGA, m=6.0)
    large = _cdfs(model, params, PGA, m=7.5)
    assert np.all(large > small)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_cdf_decreases_with_penetration_resistance(model):
    params, resistance = MODELS[model]
    stronger = copy.deepcopy(params)
    stronger[resistance] = _scale(params[resistance], 1.5)
    assert np.all(_cdfs(model, stronger, PGA) < _cdfs(model, params, PGA))


@pytest.mark.parametrize("model", sorted(MODELS))
def test_cdf_approaches_step_function_for_extreme_shaking(model):
    params, _ = MODELS[model]
    cdfs = _cdfs(model, params, np.array([0.001, 5.0]), sigma=0.01)
    assert np.all(cdfs[0] < 0.01)
    # ngl_smt_2024 is the probability of profile manifestation, which stays slightly below 1
    assert np.all(cdfs[1] > 0.95)
