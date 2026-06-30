import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.liquefaction_models import ngl_smt_2024

from constants import PGA_INPUTS


LIQ_CASES = [
    (
        "cetin_et_al_2018",
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "vs12": 180.0,
            "depth": 5.0,
            "n160": 12.0,
            "fc": 8.0,
            "pa": 101.325,
        },
    ),
    (
        "cetin_et_al_2018",
        {
            "sigmav": 300.0,
            "sigmavp": 200.0,
            "vs12": 180.0,
            "depth": 25.0,
            "n160": 12.0,
            "fc": 8.0,
            "pa": 101.325,
        },
    ),
    (
        "moss_et_al_2006",
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "depth": 5.0,
            "qc": 5.0,
            "fs": 0.5,
            "pa": 101.325,
        },
    ),
    (
        "moss_et_al_2006",
        {
            "sigmav": 300.0,
            "sigmavp": 200.0,
            "depth": 25.0,
            "qc": 20.0,
            "fs": 0.01,
            "pa": 101.325,
        },
    ),
    (
        "boulanger_idriss_2016",
        {"sigmav": 40.0, "sigmavp": 28.0, "depth": 2.2, "qc1ncs": 71.8, "pa": 101.325},
    ),
    (
        "boulanger_idriss_2012",
        {
            "sigmav": 100.0,
            "sigmavp": 60.0,
            "depth": 5.0,
            "n160": 14.0,
            "fc": 5.0,
            "pa": 101.325,
        },
    ),
    (
        "ngl_smt_2024",
        {
            "ztop": [0.0, 1.0, 5.0],
            "zbot": [1.0, 5.0, 9.0],
            "qc1ncs": [5.0, 150.0, 400.0],
            "ic": [1.5, 2.0, 2.6],
            "sigmav": [20.0, 60.0, 120.0],
            "sigmavp": [15.0, 45.0, 90.0],
            "ksat": [1.0, 1.0, 1.0],
            "pa": 101.325,
        },
    ),
]


@pytest.mark.parametrize("model,params", LIQ_CASES)
def test_get_liquefaction_cdfs(model, params):
    m, mu_pga, sigma_pga, fsl = PGA_INPUTS
    config = {"liquefaction_models": {model: params}}

    cdfs, eps = plha.get_liquefaction_cdfs(m, mu_pga, sigma_pga, fsl, model, config)

    assert cdfs.shape[-1] == len(fsl)
    assert np.all(np.isfinite(cdfs))
    assert np.all((cdfs >= 0.0) & (cdfs <= 1.0))


def test_get_liquefaction_cdfs_unknown_model_returns_none():
    m, mu_pga, sigma_pga, fsl = PGA_INPUTS
    result = plha.get_liquefaction_cdfs(m, mu_pga, sigma_pga, fsl, "not_a_model", {})
    assert result is None


def test_ngl_pdf_is_unit_normal_at_peak():
    assert ngl_smt_2024.pdf(0.0, 0.0, 1.0) == pytest.approx(1.0 / np.sqrt(2.0 * np.pi))
