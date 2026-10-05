"""
Compare ground motion models against reference values from pygmm, an independent implementation
of the NGA-West2 models (ASK14, BSSA14, CB14, CY14, and Idriss 2014). See tests/data/make_gmm_reference.py for how the reference values were
generated. Scenarios cover reverse, normal, and strike-slip ruptures, M5 to M7.8, footwall and
hanging wall sites from 2 km to 150 km, three Vs30 values, and measured and inferred Vs30.
A second set of scenarios checks the ASK14 soil depth term for user-specified z1p0.
"""

import itertools
import json
import warnings
from pathlib import Path

import numpy as np
import pytest

from ucla_plha import plha

_DATA = json.loads((Path(__file__).parent / "data" / "gmm_reference.json").read_text())
REFERENCE = _DATA["scenarios"]
ASK14_BASIN = _DATA["ask14_basin_scenarios"]
GMMS = ["bssa14", "ask14", "cb14", "cy14", "idriss14"]
# Idriss (2014) warns for Vs30 outside 450 to 1200 m/s. The reference scenarios include softer
# sites to test the model equations over the range of Vs30 used for liquefaction hazard.
IDRISS_VS30_WARNING = "ignore:vs30 = .* idriss14:UserWarning"
VS30_CASES = sorted({(s["vs30"], s["measured_vs30"]) for s in REFERENCE})
# Tolerance on natural log of PGA and on its standard deviation
ATOL = 1.0e-3


def _scenarios(vs30, measured_vs30, select=lambda s: True, reference=REFERENCE):
    rows = [
        s
        for s in reference
        if s["vs30"] == vs30 and s["measured_vs30"] == measured_vs30 and select(s)
    ]
    keys = ["fault_type", "m", "rjb", "rrup", "rx", "rx1", "ry0", "dip", "ztor", "zbor"]
    arrays = {k: np.array([s[k] for s in rows]) for k in keys}
    arrays["fault_type"] = arrays["fault_type"].astype(int)
    return rows, arrays


def _run(gmm, vs30, measured_vs30, a, z1p0=None):
    return plha.get_ground_motion_data(
        gmm,
        vs30,
        measured_vs30,
        z1p0,
        None,
        a["fault_type"],
        a["rjb"],
        a["rrup"],
        a["rx"],
        a["rx1"],
        a["ry0"],
        a["m"],
        a["ztor"],
        a["zbor"],
        a["dip"],
    )


def _describe(rows, i):
    s = rows[i]
    return (
        f"fault_type={s['fault_type']} M={s['m']} Rx={s['rx']} Rrup={s['rrup']:.2f} "
        f"dip={s['dip']} ztor={s['ztor']} Vs30={s['vs30']} measured={s['measured_vs30']}"
    )


@pytest.mark.filterwarnings(IDRISS_VS30_WARNING)
@pytest.mark.parametrize(
    "gmm,vs30,measured_vs30",
    [(g, v, mv) for g, (v, mv) in itertools.product(GMMS, VS30_CASES)],
)
def test_gmm_matches_pygmm(gmm, vs30, measured_vs30):
    rows, a = _scenarios(vs30, measured_vs30)
    mu, sigma = _run(gmm, vs30, measured_vs30, a)
    mu_ref = np.array([s["expected"][gmm]["mu"] for s in rows])
    sigma_ref = np.array([s["expected"][gmm]["sigma"] for s in rows])

    i = np.argmax(np.abs(mu - mu_ref))
    assert (
        np.abs(mu - mu_ref)[i] < ATOL
    ), f"{gmm} ln(PGA) = {mu[i]:.4f}, pygmm = {mu_ref[i]:.4f} for {_describe(rows, i)}"
    i = np.argmax(np.abs(sigma - sigma_ref))
    assert (
        np.abs(sigma - sigma_ref)[i] < ATOL
    ), f"{gmm} sigma = {sigma[i]:.4f}, pygmm = {sigma_ref[i]:.4f} for {_describe(rows, i)}"


@pytest.mark.parametrize("vs30", sorted({s["vs30"] for s in ASK14_BASIN}))
def test_ask14_basin_term_matches_pygmm(vs30):
    for z1p0 in sorted({s["z1p0"] for s in ASK14_BASIN if s["vs30"] == vs30}):
        rows, a = _scenarios(
            vs30, False, lambda s: s["z1p0"] == z1p0, reference=ASK14_BASIN
        )
        mu, sigma = _run("ask14", vs30, False, a, z1p0=z1p0)
        mu_ref = np.array([s["expected"]["ask14"]["mu"] for s in rows])
        sigma_ref = np.array([s["expected"]["ask14"]["sigma"] for s in rows])
        i = np.argmax(np.abs(mu - mu_ref))
        assert np.abs(mu - mu_ref)[i] < ATOL, (
            f"ask14 ln(PGA) = {mu[i]:.4f}, pygmm = {mu_ref[i]:.4f} for z1p0={z1p0} km, "
            f"{_describe(rows, i)}"
        )
        np.testing.assert_allclose(sigma, sigma_ref, atol=ATOL)


@pytest.mark.parametrize("vs30", [180.0, 200.0, 270.0, 400.0, 600.0, 900.0])
def test_ask14_default_z1p0_matches_unspecified(vs30):
    # Specifying z1p0 equal to the default Z1ref(Vs30) from ASK14 Equation 18 must give the same
    # result as leaving z1p0 unspecified
    _, a = _scenarios(250.0, False)
    z1ref = np.exp(-7.67 / 4 * np.log((vs30**4 + 610**4) / (1360**4 + 610**4))) / 1000
    mu, sigma = _run("ask14", vs30, False, a)
    mu_z1, sigma_z1 = _run("ask14", vs30, False, a, z1p0=z1ref)
    np.testing.assert_allclose(mu_z1, mu, atol=1e-12)
    np.testing.assert_allclose(sigma_z1, sigma, atol=1e-12)


def test_ask14_basin_term_is_continuous_in_vs30():
    # The soil depth coefficients are interpolated between bin centers, so PGA should change
    # smoothly with Vs30 across the bin edges at 200, 300, and 500 m/s for a deep basin site
    _, a = _scenarios(250.0, False)
    for edge in (200.0, 300.0, 500.0):
        below, _ = _run("ask14", edge - 0.01, False, a, z1p0=1.5)
        above, _ = _run("ask14", edge + 0.01, False, a, z1p0=1.5)
        np.testing.assert_allclose(above, below, atol=1e-3)


@pytest.mark.parametrize("gmm", ["ask14", "cb14", "cy14"])
def test_hanging_wall_amplifies_dipping_ruptures(gmm):
    # For a large dipping rupture, a site above the rupture on the hanging wall should have higher
    # PGA than a footwall site at the same Rrup. This isolates the hanging wall term.
    m = np.array([7.0, 7.0])
    dip = np.array([45.0, 45.0])
    ztor = np.array([0.0, 0.0])
    zbor = np.array([15.0, 15.0])
    rrup = np.array([7.07, 7.07])
    rjb = np.array([0.0, 7.07])
    rx = np.array([10.0, -7.07])
    rx1 = rx - 15.0
    mu, _ = plha.get_ground_motion_data(
        gmm,
        760.0,
        False,
        None,
        None,
        np.array([1, 1]),
        rjb,
        rrup,
        rx,
        rx1,
        np.zeros(2),
        m,
        ztor,
        zbor,
        dip,
    )
    assert mu[0] > mu[1] + 0.1


@pytest.mark.filterwarnings(IDRISS_VS30_WARNING)
@pytest.mark.parametrize("gmm", GMMS)
def test_pga_decreases_with_distance_and_increases_with_magnitude(gmm):
    n = 20
    r = np.logspace(0, np.log10(200), n)
    base = {
        "fault_type": np.full(n, 3),
        "rjb": r,
        "rrup": np.hypot(r, 2.0),
        "rx": r,
        "rx1": r,
        "ry0": np.zeros(n),
        "dip": np.full(n, 90.0),
        "ztor": np.full(n, 2.0),
        "zbor": np.full(n, 14.0),
    }
    mu = {}
    for m in (5.5, 6.5, 7.5):
        mu[m], _ = _run(gmm, 400.0, False, dict(base, m=np.full(n, m)))
        assert np.all(np.diff(mu[m]) < 0.0)
    assert np.all(mu[6.5] > mu[5.5]) and np.all(mu[7.5] > mu[6.5])


@pytest.mark.parametrize(
    "vs30,warns", [(300.0, True), (450.0, False), (1200.0, False), (1500.0, True)]
)
def test_idriss14_warns_outside_recommended_vs30(vs30, warns):
    _, a = _scenarios(250.0, False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _run("idriss14", vs30, False, a)
    messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
    assert any("idriss14" in msg for msg in messages) == warns


@pytest.mark.filterwarnings(IDRISS_VS30_WARNING)
def test_idriss14_treats_normal_faulting_as_strike_slip():
    # Idriss (2014) has a style of faulting term for reverse earthquakes only
    _, a = _scenarios(760.0, False)
    a_ss = dict(a, fault_type=np.full(len(a["m"]), 3))
    a_ns = dict(a, fault_type=np.full(len(a["m"]), 2))
    a_rs = dict(a, fault_type=np.full(len(a["m"]), 1))
    mu_ss, _ = _run("idriss14", 760.0, False, a_ss)
    mu_ns, _ = _run("idriss14", 760.0, False, a_ns)
    mu_rs, _ = _run("idriss14", 760.0, False, a_rs)
    np.testing.assert_array_equal(mu_ns, mu_ss)
    np.testing.assert_allclose(mu_rs - mu_ss, 0.08, atol=1e-12)
