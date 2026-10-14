"""Input parameter uncertainty of Cetin et al. (2018b) in the cetin_et_al_2018 model."""

import json
import math
import warnings
from importlib.resources import files

import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.liquefaction_models import cetin_et_al_2018 as cetin

SITE = dict(sigmav=100.0, sigmavp=60.0, vs12=180.0, depth=5.0, n160=12.0, fc=8.0, pa=101.325)
UNCERTAINTIES = dict(sigma_n160=2.3, sigma_fc=5.3, cov_sigmavp=0.1, cov_csr=0.1)
M = np.array([5.5, 6.5, 7.5])


def test_default_is_model_uncertainty_only():
    _, sigma = cetin.get_ln_crr(M, SITE["n160"], SITE["fc"], SITE["sigmavp"], SITE["pa"])
    np.testing.assert_array_equal(sigma, 2.95 / 11.771)


def test_hand_calculation():
    n160, fc = SITE["n160"], SITE["fc"]
    # Eqs. 4 and 5 of Cetin et al. (2018b), written out
    sigma_input_sq = (
        11.771**2 * 0.1**2
        + 2.3**2 * (1 + 0.00167 * 8.0) ** 2
        + 5.3**2 * (0.00167 * 12.0 + 0.089) ** 2
        + 3.958**2 * 0.1**2
    )
    expected = math.sqrt(0.392 * sigma_input_sq + 2.95**2)
    assert cetin.get_sigma_total(n160, fc, True, **UNCERTAINTIES) == pytest.approx(expected, rel=1e-12)
    assert expected == pytest.approx(3.4011, abs=1e-4)
    mu0, sigma0 = cetin.get_ln_crr(M, n160, fc, SITE["sigmavp"], SITE["pa"])
    mu1, sigma1 = cetin.get_ln_crr(M, n160, fc, SITE["sigmavp"], SITE["pa"], True, **UNCERTAINTIES)
    # The parameter uncertainty only changes the standard deviation, not the mean
    np.testing.assert_array_equal(mu1, mu0)
    np.testing.assert_allclose(sigma1, expected / 11.771, rtol=1e-12)


def test_zero_uncertainties_equal_the_default():
    assert cetin.get_sigma_total(12.0, 8.0, True) == 2.95


def _cdfs(c, fsl=np.array([0.8, 1.0, 1.2])):
    config = {"liquefaction_models": {"cetin_et_al_2018": dict(SITE, weight=1.0, **c)}}
    mu_ln_pga = np.log(np.array([0.1, 0.2, 0.4]))
    sigma_ln_pga = np.full(3, 0.6)
    return plha.get_liquefaction_cdfs(M, mu_ln_pga, sigma_ln_pga, fsl, "cetin_et_al_2018", config)


def test_config_options():
    default, _ = _cdfs({})
    off, _ = _cdfs({"parameter_uncertainty": False, **UNCERTAINTIES})
    np.testing.assert_array_equal(off, default)
    on, _ = _cdfs({"parameter_uncertainty": True, **UNCERTAINTIES})
    assert np.all(on != default)
    # A larger standard deviation moves the probabilities toward 0.5
    assert np.all(np.abs(on - 0.5) <= np.abs(default - 0.5))


def test_warning_when_enabled_without_uncertainties():
    default, _ = _cdfs({})
    with pytest.warns(UserWarning, match="parameter_uncertainty"):
        on, _ = _cdfs({"parameter_uncertainty": True})
    np.testing.assert_array_equal(on, default)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _cdfs({"parameter_uncertainty": True, "sigma_n160": 2.0})


def test_schema_has_the_options():
    schema = json.loads(files("ucla_plha").joinpath("ucla_plha_schema.json").read_text())
    props = schema["properties"]["liquefaction_models"]["properties"]["cetin_et_al_2018"]["properties"]
    assert props["parameter_uncertainty"]["type"] == "boolean"
    for key in UNCERTAINTIES:
        assert props[key]["type"] == "number"
