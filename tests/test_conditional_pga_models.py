"""Boulanger and Idriss (2016) and Moss et al. (2006) conditioned on PGA and integrated on the
magnitude x ln(PGA) rate grid, and their earlier closed forms ("pga_uncertainty": "closed_form")."""

import copy

import jsonschema
import numpy as np
import pytest
from scipy.special import ndtr
from scipy.stats import norm

from ucla_plha import plha
from ucla_plha.liquefaction_models import boulanger_idriss_2016, moss_et_al_2006
from pga_grid_helpers import BINS, CONFIG, run_hazard, validate

FSL = np.array([0.3, 0.5, 1.0, 2.0])
MODELS = {
    "boulanger_idriss_2016": {"depth": 2.2, "qc1ncs": 71.8, "sigmav": 40.0, "sigmavp": 28.0, "pa": 101.325},
    # Rf = 0.48%: no delta qc adjustment
    "moss_et_al_2006": {"depth": 2.85, "qc": 5.49, "fs": 0.02617, "sigmav": 46.55, "sigmavp": 43.22, "pa": 101.325},
}
# Moss et al. (2006) with Rf = 2% (delta qc = x1 ln CSR + x2)
MOSS_FINES = {"depth": 5.0, "qc": 5.0, "fs": 0.1, "sigmav": 100.0, "sigmavp": 60.0, "pa": 101.325}
MODULES = {"boulanger_idriss_2016": boulanger_idriss_2016, "moss_et_al_2006": moss_et_al_2006}


def _args(model, params):
    if model == "boulanger_idriss_2016":
        keys = ("sigmav", "sigmavp", "depth", "qc1ncs")
    else:
        keys = ("sigmav", "sigmavp", "depth", "qc", "fs")
    return [params[k] for k in keys]


def _events(n=300, seed=3):
    rng = np.random.default_rng(seed)
    return rng.uniform(5.0, 8.0, n), np.log(rng.uniform(0.01, 1.0, n)), rng.uniform(0.4, 0.8, n)


def _old_bi16(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, qc1ncs, fsl, pa):
    """boulanger_idriss_2016.get_fsl_cdfs of ucla_plha 2.x."""
    mu_ln_crr, sigma_ln_crr = boulanger_idriss_2016.get_ln_crr(m, qc1ncs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = boulanger_idriss_2016.get_ln_csr(
        mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth
    )
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    std_ln_fsl = np.sqrt(sigma_ln_csr**2 + sigma_ln_crr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / std_ln_fsl[:, np.newaxis]
    return ndtr(eps), eps


def _old_moss(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, qc, fs, fsl, pa):
    """moss_et_al_2006.get_fsl_cdfs closed form of this version: ucla_plha 2.x without the duration
    weighting (CSR = 0.65 amax (sigmav / sigmav') rd at the median PGA)."""
    rd = moss_et_al_2006.get_rd(mu_ln_pga, m, depth)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    mu_ln_crr, sigma_ln_crr = moss_et_al_2006.get_ln_crr(m, mu_ln_csr, qc, fs, sigmavp, pa)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_pga**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]
    return ndtr(eps), eps


def _moss_2x(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, qc, fs, fsl, pa):
    """moss_et_al_2006.get_fsl_cdfs of ucla_plha 2.x (CSR / DWF), for comparison only."""
    dwf = 17.84 * np.clip(m, 5.5, 8.5) ** -1.43
    rd = moss_et_al_2006.get_rd(mu_ln_pga, m, depth)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd) - np.log(dwf)
    mu_ln_crr, sigma_ln_crr = moss_et_al_2006.get_ln_crr(m, mu_ln_csr, qc, fs, sigmavp, pa)
    eps = (np.log(fsl) - (mu_ln_crr - mu_ln_csr)[:, np.newaxis]) / np.sqrt(
        sigma_ln_crr**2 + sigma_ln_pga**2
    )[:, np.newaxis]
    return ndtr(eps), eps


OLD = {"boulanger_idriss_2016": _old_bi16, "moss_et_al_2006": _old_moss}
CASES = [(model, params) for model, params in MODELS.items()] + [("moss_et_al_2006", MOSS_FINES)]


@pytest.mark.parametrize("model,params", CASES)
def test_closed_form_reproduces_frozen_cdfs_bit_for_bit(model, params):
    # BI16: ucla_plha 2.x; Moss et al. (2006): 2.x without the DWF (_old_moss)
    m, mu, sigma = _events()
    args = (mu, sigma, m, *_args(model, params), FSL, params["pa"])
    cdfs, eps = MODULES[model].get_fsl_cdfs(*args, pga_uncertainty="closed_form")
    old_cdfs, old_eps = OLD[model](*args)
    assert np.array_equal(cdfs, old_cdfs)
    assert np.array_equal(eps, old_eps)


@pytest.mark.parametrize("truncation_level", [None, 3.0])
@pytest.mark.parametrize("model,params", CASES)
def test_conditional_quadrature_matches_brute_force(model, params, truncation_level):
    m, mu, sigma = np.array([6.8]), np.array([np.log(0.25)]), np.array([0.6])
    z_hi = 10.0 if truncation_level is None else truncation_level
    z = np.linspace(-10.0, z_hi, 200001)
    p = MODULES[model].get_conditional_fsl_cdfs(
        np.full(len(z), 6.8), mu[0] + sigma[0] * z, *_args(model, params), FSL, params["pa"]
    )
    weight = norm.pdf(z) / ndtr(z_hi)
    expected = np.trapezoid(p * weight[:, None], z, axis=0)
    cdfs, eps = MODULES[model].get_fsl_cdfs(
        mu, sigma, m, *_args(model, params), FSL, params["pa"],
        truncation_level=truncation_level, n_quadrature=256,
    )
    # the rd of Moss et al. (2006) has a kink at the amax cap (1 g), which limits the
    # Gauss-Legendre accuracy
    rtol = 1e-6
    if model == "moss_et_al_2006":
        rtol = 5e-4 if truncation_level is None else 2e-5
    np.testing.assert_allclose(cdfs[0], expected, rtol=rtol, atol=1e-12)
    np.testing.assert_allclose(eps, norm.ppf(cdfs))


def test_bi16_conditional_untruncated_equals_closed_form():
    # PGA enters only through CSR, so the integral over the untruncated distribution is the
    # closed form
    params = MODELS["boulanger_idriss_2016"]
    m, mu, sigma = _events()
    args = (mu, sigma, m, *_args("boulanger_idriss_2016", params), FSL, params["pa"])
    closed, _ = boulanger_idriss_2016.get_fsl_cdfs(*args, pga_uncertainty="closed_form")
    conditional, _ = boulanger_idriss_2016.get_fsl_cdfs(*args)
    np.testing.assert_allclose(conditional, closed, rtol=1e-9, atol=1e-10)


def test_moss_conditional_defined_for_small_csr():
    # with Rf > 0.5, delta qc = x1 ln CSR + x2 makes qc1 + delta qc negative for very small
    # CSR (nan in the 2.x formula); the conditional probability clips it at 0
    params = dict(MOSS_FINES, fs=0.5)  # Rf = 10% (clipped to 5%)
    x = np.log(np.array([1e-6, 1e-4, 0.01, 0.3]))
    p = moss_et_al_2006.get_conditional_fsl_cdfs(
        np.full(4, 7.0), x, *_args("moss_et_al_2006", params), FSL, params["pa"]
    )
    assert np.all(np.isfinite(p))
    assert np.all(p[:2] < 1e-12) and np.all(np.diff(p, axis=0) >= 0.0)


def test_moss_delta_qc_nonnegative():
    # Rf = 2%: delta qc = x1 ln CSR + x2 >= 0 in the conditional path; x1 ln CSR + x2 < 0 only for
    # CSR < exp(-x2 / x1) = exp(-3.84)
    p = MOSS_FINES
    rf = p["fs"] / p["qc"] * 100
    x1, x2 = 0.38 * rf - 0.19, 1.46 * rf - 0.73
    assert np.isclose(-x2 / x1, -3.8421, atol=1e-4)
    ln_csr = np.array([-8.0, -5.0, -x2 / x1, -3.0, -1.0])
    m = np.full(len(ln_csr), 7.0)
    args = (p["qc"], p["fs"], p["sigmavp"], p["pa"])
    clipped, _ = moss_et_al_2006.get_ln_crr(m, ln_csr, *args, nonnegative_delta_qc=True)
    formula, _ = moss_et_al_2006.get_ln_crr(m, ln_csr, *args)
    at_zero, _ = moss_et_al_2006.get_ln_crr(m, np.full(len(ln_csr), -x2 / x1), *args)
    np.testing.assert_allclose(clipped[:3], at_zero[:3], rtol=1e-12)
    assert np.all(formula[:2] < clipped[:2])
    assert np.array_equal(clipped[3:], formula[3:])
    # the hazard with Rf > 0.5% is finite and below the closed form at small CSR
    curve = _plha(run_hazard(_config("moss_et_al_2006", MOSS_FINES)))
    assert np.all(np.isfinite(curve)) and np.all(curve > 0.0)


def test_moss_no_duration_weighting():
    # CSR = 0.65 amax (sigmav / sigmav') rd, without the DWF of ucla_plha 2.x, in every treatment
    params = MOSS_FINES
    x = np.log(np.array([0.05, 0.2, 0.6]))
    m = np.array([6.0, 7.0, 8.0])
    csr, _ = moss_et_al_2006.get_ln_csr(x, 0.0, m, params["sigmav"], params["sigmavp"], params["depth"])
    rd = moss_et_al_2006.get_rd(x, m, params["depth"])
    np.testing.assert_allclose(csr, x + np.log(0.65 * params["sigmav"] / params["sigmavp"] * rd), rtol=1e-14)
    assert not hasattr(moss_et_al_2006, "DURATION_WEIGHTING")
    # the closed form differs from 2.x (CSR / DWF), which can't be reproduced
    mm, mu, sigma = _events()
    args = (mu, sigma, mm, *_args("moss_et_al_2006", params), FSL, params["pa"])
    closed, _ = moss_et_al_2006.get_fsl_cdfs(*args, pga_uncertainty="closed_form")
    assert not np.allclose(closed, _moss_2x(*args)[0], rtol=1e-3)


@pytest.mark.parametrize("treatment", ["conditional", "closed_form"])
def test_moss_duration_weighting_rejected(treatment):
    config = _config("moss_et_al_2006", pga_uncertainty=treatment, duration_weighting=False)
    with pytest.raises(ValueError, match="duration_weighting"):
        run_hazard(config)
    with pytest.raises(ValueError, match="double counted"):
        plha.get_liquefaction_cdfs(
            np.array([7.0]), np.array([-1.0]), np.array([0.6]), FSL, "moss_et_al_2006",
            {"liquefaction_models": {"moss_et_al_2006": dict(MODELS["moss_et_al_2006"], duration_weighting=True)}},
        )


def test_moss_rd_amax_cap():
    # rd of the conditional probability uses amax at most 1 g (RD_AMAX_CAP); without the cap, rd
    # drops to 0.2 and jumps to 1.0 between about 2.5 and 4 g, so P(FS < fsl | PGA) decreases
    params = MODELS["moss_et_al_2006"]
    x = np.log(np.geomspace(0.05, 20.0, 60))
    m = np.full(len(x), 7.0)
    rd = moss_et_al_2006.get_rd(x, m, params["depth"], 1.0)
    assert np.array_equal(rd, moss_et_al_2006.get_rd(np.minimum(x, 0.0), m, params["depth"]))
    assert moss_et_al_2006.RD_AMAX_CAP == 1.0
    capped = moss_et_al_2006.get_conditional_fsl_cdfs(m, x, *_args("moss_et_al_2006", params), FSL, params["pa"])
    uncapped = moss_et_al_2006.get_conditional_fsl_cdfs(
        m, x, *_args("moss_et_al_2006", params), FSL, params["pa"], rd_amax_cap=None
    )
    assert np.all(np.diff(capped, axis=0) >= 0.0)
    assert np.any(np.diff(uncapped, axis=0) < -0.1)
    np.testing.assert_array_equal(capped[np.exp(x) <= 1.0], uncapped[np.exp(x) <= 1.0])


def test_moss_rd_amax_cap_config():
    default = _plha(run_hazard(_config("moss_et_al_2006")))
    assert np.array_equal(default, _plha(run_hazard(_config("moss_et_al_2006", rd_amax_cap=1.0))))
    uncapped = _plha(run_hazard(_config("moss_et_al_2006", rd_amax_cap=None)))
    assert not np.array_equal(default, uncapped)
    # closed_form (rd of the median PGA) has no cap
    closed = _plha(run_hazard(_config("moss_et_al_2006", pga_uncertainty="closed_form")))
    closed_cap = _plha(run_hazard(_config("moss_et_al_2006", pga_uncertainty="closed_form", rd_amax_cap=0.1)))
    assert np.array_equal(closed, closed_cap)


def _config(model, params=None, truncation_level=None, **opts):
    config = copy.deepcopy(CONFIG)
    config["output"]["plha"]["fsl"] = FSL.tolist()
    config["liquefaction_models"] = {model: dict(params or MODELS[model], weight=1.0, **opts)}
    if truncation_level is not None:
        config["constraints"]["truncation_level"] = truncation_level
    return config


def _plha(out):
    return np.array(out["output"]["plha"]["annual_rate_of_nonexceedance"])


@pytest.mark.parametrize("truncation_level", [None, 3.0])
@pytest.mark.parametrize("model,params", CASES)
def test_hazard_grid_matches_quadrature_and_disaggregation(model, params, truncation_level):
    # Moss et al. (2006): the kink of rd at the amax cap (1 g) limits both the grid and the
    # 128-point quadrature to about 1e-3 at FS_L 0.3 (both checked against a fine integration)
    rtol = 2e-3 if model == "moss_et_al_2006" else 5e-4
    grid = _plha(run_hazard(_config(model, params, truncation_level)))
    quad = _plha(run_hazard(_config(model, params, truncation_level, pga_integration="quadrature")))
    np.testing.assert_allclose(grid, quad, rtol=rtol)
    config = _config(model, params, truncation_level)
    config["output"]["plha"]["disaggregation"] = BINS
    out = run_hazard(config)
    np.testing.assert_allclose(_plha(out), quad, rtol=rtol)
    disagg = np.array(out["output"]["plha"]["disaggregation"])
    np.testing.assert_allclose(disagg.reshape(len(FSL), -1).sum(axis=1), 100.0, rtol=1e-10)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_truncation(model):
    curves = {}
    for treatment in ("conditional", "closed_form"):
        for truncation_level in (None, 2.0):
            config = _config(model, truncation_level=truncation_level, pga_uncertainty=treatment)
            curves[treatment, truncation_level] = _plha(run_hazard(config))
    assert np.array_equal(curves["closed_form", None], curves["closed_form", 2.0])
    assert np.all(curves["conditional", 2.0] < curves["conditional", None])


@pytest.mark.parametrize("model", sorted(MODELS))
def test_hazard_closed_form_reproduces_old_model(model, monkeypatch):
    new = _plha(run_hazard(_config(model, pga_uncertainty="closed_form")))

    def old(*args, **kwargs):
        return OLD[model](*args)

    monkeypatch.setattr(MODULES[model], "get_fsl_cdfs", old)
    assert np.array_equal(new, _plha(run_hazard(_config(model, pga_uncertainty="closed_form"))))


def test_default_switch(monkeypatch):
    # DEFAULT_PGA_UNCERTAINTY gives the treatment of a model without "pga_uncertainty"
    closed = _plha(run_hazard(_config("moss_et_al_2006", pga_uncertainty="closed_form")))
    monkeypatch.setitem(plha.DEFAULT_PGA_UNCERTAINTY, "moss_et_al_2006", "closed_form")
    assert np.array_equal(_plha(run_hazard(_config("moss_et_al_2006"))), closed)


def test_schema_options():
    for model in MODELS:
        good = _config(model, pga_uncertainty="closed_form", pga_integration="quadrature", pga_grid={"d_m": 0.1})
        validate(good)
        with pytest.raises(jsonschema.ValidationError):
            validate(_config(model, pga_uncertainty="per_layer"))

