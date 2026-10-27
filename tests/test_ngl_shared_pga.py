"""NGL SMT 2024 with the peak acceleration shared by the layers of a profile, and the magnitude x
ln(PGA) rate grid (ucla_plha.pga_grid)."""

import copy

import jsonschema
import numpy as np
import pytest
from scipy.special import ndtr
from scipy.stats import norm

from ucla_plha import pga_grid, plha
from ucla_plha.liquefaction_models import ngl_smt_2024
from pga_grid_helpers import BINS, CONFIG, NGL, run_hazard, validate

FSL = np.array([0.5, 1.0, 2.0])


def _events(n=300, seed=3):
    rng = np.random.default_rng(seed)
    m = rng.uniform(5.0, 8.0, n)
    mu = np.log(rng.uniform(0.005, 1.0, n))
    sigma = rng.uniform(0.4, 0.8, n)
    return m, mu, sigma


def _old_get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa):
    """ngl_smt_2024.get_fsl_cdfs of ucla_plha 2.x (before the "pga_uncertainty" option)."""
    tc = 2.0
    mu_ln_crr, sigma_ln_crr = ngl_smt_2024.get_ln_crr(mu_ln_pga, m, qc1Ncs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = ngl_smt_2024.get_ln_csr(
        mu_ln_pga, sigma_ln_pga, m, ztop, zbot, sigmav, sigmavp
    )
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps_triggering = (
        np.log(fsl[np.newaxis, np.newaxis, :]) - mu_ln_fsl[:, :, np.newaxis]
    ) / sigma_ln_fsl[:, np.newaxis, np.newaxis]
    pfts = ndtr(eps_triggering)
    pfs = 1.0 - 1.0 / (1.0 + np.exp(-14.8 * (Ic / 2.635 - 1.0)))
    pfmt = 1.0 / (1.0 + np.exp(-(7.613 - 0.338 * ztop - 3.042 * Ic)))
    t = zbot - ztop
    fsl_cdfs = 1.0 - np.prod(
        (
            1.0
            - pfmt[np.newaxis, :, np.newaxis]
            * pfts
            * pfs[np.newaxis, :, np.newaxis]
            * Ksat[np.newaxis, :, np.newaxis]
        )
        ** (t[np.newaxis, :, np.newaxis] / tc),
        axis=1,
    )
    return fsl_cdfs, norm.ppf(fsl_cdfs)


def test_per_layer_reproduces_old_cdfs_bit_for_bit():
    m, mu, sigma = _events()
    args = (mu, sigma, m, *NGL.values(), FSL, 101.325)
    cdfs, eps = ngl_smt_2024.get_fsl_cdfs(*args, pga_uncertainty="per_layer")
    old_cdfs, old_eps = _old_get_fsl_cdfs(*args)
    assert np.array_equal(cdfs, old_cdfs)
    assert np.array_equal(eps, old_eps)


def test_shared_quadrature_matches_brute_force_integration():
    # one event, the conditional probability integrated with a fine trapezoidal rule
    m, mu, sigma = np.array([7.0]), np.array([np.log(0.3)]), np.array([0.6])
    z = np.linspace(-10.0, 10.0, 40001)
    p = ngl_smt_2024.get_conditional_fsl_cdfs(
        np.full(len(z), 7.0), mu[0] + sigma[0] * z, *NGL.values(), FSL, 101.325
    )
    expected = np.trapezoid(p * norm.pdf(z)[:, None], z, axis=0)
    cdfs, _ = ngl_smt_2024.get_fsl_cdfs(mu, sigma, m, *NGL.values(), FSL, 101.325)
    np.testing.assert_allclose(cdfs[0], expected, rtol=1e-7, atol=1e-12)
    # the shared PGA gives lower probabilities than independent PGAs in the layers
    old, _ = ngl_smt_2024.get_fsl_cdfs(
        mu, sigma, m, *NGL.values(), FSL, 101.325, pga_uncertainty="per_layer"
    )
    assert np.all(cdfs <= old)


def test_conditional_cdfs_single_layer_closed_form():
    # one saturated layer: P = 1 - (1 - pfmt pfs Phi(...))^(t / tc)
    layer = {k: v[1:2] for k, v in NGL.items()}
    m, x = np.array([6.5, 7.5]), np.log(np.array([0.2, 0.4]))
    p = ngl_smt_2024.get_conditional_fsl_cdfs(m, x, *layer.values(), FSL, 101.325)
    mu_crr, _ = ngl_smt_2024.get_ln_crr(x, m, layer["qc1Ncs"], layer["sigmavp"], 101.325)
    mu_csr, _ = ngl_smt_2024.get_ln_csr(x, 0.0, m, layer["ztop"], layer["zbot"], layer["sigmav"], layer["sigmavp"])
    pfs = 1.0 - 1.0 / (1.0 + np.exp(-14.8 * (layer["Ic"] / 2.635 - 1.0)))
    pfmt = 1.0 / (1.0 + np.exp(-(7.613 - 0.338 * layer["ztop"] - 3.042 * layer["Ic"])))
    pfts = ndtr((np.log(FSL)[None, :] - (mu_crr - mu_csr)) / 0.44)
    t = layer["zbot"] - layer["ztop"]
    expected = 1.0 - (1.0 - pfmt * pfs * pfts) ** (t / 2.0)
    np.testing.assert_allclose(p, expected, rtol=1e-13)


@pytest.mark.parametrize("truncation_level", [None, 3.0, 1.5])
def test_bin_weights_exact_for_quadratics(truncation_level):
    _, mu, sigma = _events(50)
    mu = np.log(np.random.default_rng(4).uniform(0.02, 0.3, 50))
    # (the first and last bins hold the tails: limits far from the means)
    grid = pga_grid.RateGrid(truncation_level, pga_min=1e-6, pga_max=1e3)
    W = grid.event_weights(mu, sigma)
    x = grid.ln_pga
    if truncation_level is None:
        mean, second = mu, mu**2 + sigma**2
    else:
        t = truncation_level
        lam = norm.pdf(t) / ndtr(t)
        mean = mu - sigma * lam
        second = mean**2 + sigma**2 * (1.0 - t * lam - lam**2)
    np.testing.assert_allclose(W.sum(axis=1), 1.0, rtol=1e-12)
    np.testing.assert_allclose(W @ x, mean, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(W @ x**2, second, rtol=1e-9)


@pytest.mark.parametrize("truncation_level", [None, 3.0])
def test_grid_and_event_cdfs_match_quadrature(truncation_level):
    m, mu, sigma = _events(2000)
    rate = np.random.default_rng(5).uniform(1e-5, 1e-3, len(m))
    args = (*NGL.values(), FSL, 101.325)
    exact, _ = ngl_smt_2024.get_fsl_cdfs(mu, sigma, m, *args, truncation_level=truncation_level, return_eps=False)
    grid = pga_grid.RateGrid(truncation_level)
    table = pga_grid.ConditionalTable(
        grid, lambda mm, xx: ngl_smt_2024.get_conditional_fsl_cdfs(mm, xx, *args)
    )
    events = table.event_cdfs(m, mu, sigma)
    np.testing.assert_allclose(events, exact, rtol=2e-3, atol=2e-5)
    np.testing.assert_allclose(rate @ events, rate @ exact, rtol=1e-4)
    for start in range(0, len(m), 700):  # in several parts
        sl = slice(start, start + 700)
        grid.add(m[sl], mu[sl], sigma[sl], rate[sl])
    np.testing.assert_allclose(table.curve(), rate @ exact, rtol=2e-4)


def test_truncation_consistent_with_psha():
    # an indicator-like conditional probability of a ground motion level reproduces the PSHA
    # exceedance (get_exceedance) with the same truncation, for a fine grid
    _, mu, sigma = _events(20)
    level = 0.2
    for truncation_level in (None, 3.0):
        grid = pga_grid.RateGrid(truncation_level, d_ln_pga=0.002, pga_min=1e-4, pga_max=20.0)
        W = grid.event_weights(mu, sigma)
        above = W[:, grid.ln_pga > np.log(level)].sum(axis=1)
        _, p = plha.get_exceedance(np.array([level]), mu, sigma, truncation_level)
        np.testing.assert_allclose(above, p[0], atol=2e-3)


def _config(**ngl):
    config = copy.deepcopy(CONFIG)
    config["liquefaction_models"]["ngl_smt_2024"].update(ngl)
    return config


def _plha(out):
    return np.array(out["output"]["plha"]["annual_rate_of_nonexceedance"])


@pytest.mark.parametrize("truncation_level", [None, 3.0])
def test_hazard_grid_matches_quadrature_and_disaggregation(truncation_level):
    config = _config()
    if truncation_level is not None:
        config["constraints"]["truncation_level"] = truncation_level
    grid = _plha(run_hazard(copy.deepcopy(config)))
    quad_config = copy.deepcopy(config)
    quad_config["liquefaction_models"]["ngl_smt_2024"]["pga_integration"] = "quadrature"
    quad = _plha(run_hazard(quad_config))
    np.testing.assert_allclose(grid, quad, rtol=1e-4)
    # disaggregation: per-event probabilities from the grid nodes
    dis_config = copy.deepcopy(config)
    dis_config["output"]["plha"]["disaggregation"] = BINS
    out = run_hazard(dis_config)
    np.testing.assert_allclose(_plha(out), quad, rtol=1e-4)
    disagg = np.array(out["output"]["plha"]["disaggregation"])
    np.testing.assert_allclose(disagg.reshape(len(FSL), -1).sum(axis=1), 100.0, rtol=1e-10)


def test_truncation_changes_shared_but_not_per_layer():
    curves = {}
    for treatment in ("shared", "per_layer"):
        for truncation_level in (None, 2.0):
            config = _config(pga_uncertainty=treatment)
            if truncation_level is not None:
                config["constraints"]["truncation_level"] = truncation_level
            curves[treatment, truncation_level] = _plha(run_hazard(config))
    assert np.array_equal(curves["per_layer", None], curves["per_layer", 2.0])
    assert np.all(curves["shared", 2.0] < curves["shared", None])


def test_hazard_per_layer_reproduces_old_model(monkeypatch):
    # get_hazard with "per_layer" gives the hazard of the ucla_plha 2.x model bit for bit
    new = _plha(run_hazard(_config(pga_uncertainty="per_layer")))

    def old(*args, return_eps=True, **kwargs):
        cdfs, eps = _old_get_fsl_cdfs(*args)
        return cdfs, (eps if return_eps else None)

    monkeypatch.setattr(ngl_smt_2024, "get_fsl_cdfs", old)
    assert np.array_equal(new, _plha(run_hazard(_config(pga_uncertainty="per_layer"))))


def test_schema_options():
    good = _config(pga_uncertainty="per_layer", pga_integration="quadrature", pga_grid={"d_ln_pga": 0.05})
    validate(good)
    with pytest.raises(jsonschema.ValidationError):
        validate(_config(pga_uncertainty="independent"))
    # grid parameters are checked by pga_grid.grid_parameters
    for pga_grid_params in ({"d_ln_pga": 0.0}, {"spacing": 0.1}, {"pga_min": 1.0, "pga_max": 0.5}):
        with pytest.raises(ValueError):
            run_hazard(_config(pga_grid=pga_grid_params))


def test_ngl_smt_2024_without_eps_same_cdfs():
    m, mu, sigma = _events(200)
    args = (mu, sigma, m, *NGL.values(), FSL, 101.325)
    for treatment in ("shared", "per_layer"):
        cdfs, eps = ngl_smt_2024.get_fsl_cdfs(*args, pga_uncertainty=treatment)
        cdfs_only, none = ngl_smt_2024.get_fsl_cdfs(*args, pga_uncertainty=treatment, return_eps=False)
        assert none is None and eps is not None
        assert np.array_equal(cdfs, cdfs_only)


def test_shared_event_cdfs_same_as_each_table():
    # the quadrature weights of the events are shared by the tables of a grid; the results are
    # those of ConditionalTable.event_cdfs of each table
    m, mu, sigma = _events(45000)  # more than one chunk
    grid = pga_grid.RateGrid(3.0)
    other = pga_grid.RateGrid(3.0, d_ln_pga=0.2)
    args = (*NGL.values(), FSL, 101.325)
    loose = dict(NGL, qc1Ncs=NGL["qc1Ncs"] * 0.8)
    tables = [
        pga_grid.ConditionalTable(grid, lambda mm, xx: ngl_smt_2024.get_conditional_fsl_cdfs(mm, xx, *args)),
        pga_grid.ConditionalTable(
            grid, lambda mm, xx: ngl_smt_2024.get_conditional_fsl_cdfs(mm, xx, *loose.values(), FSL, 101.325)
        ),
        pga_grid.ConditionalTable(other, lambda mm, xx: ngl_smt_2024.get_conditional_fsl_cdfs(mm, xx, *args)),
    ]
    shared = pga_grid.shared_event_cdfs(tables, m, mu, sigma)
    for table, p in zip(tables, shared):
        assert np.array_equal(p, table.event_cdfs(m, mu, sigma))
