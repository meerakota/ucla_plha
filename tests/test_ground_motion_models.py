import numpy as np
import pytest

from ucla_plha import plha
from ucla_plha.ground_motion_models import ask14

from constants import GMM_INPUTS


def _call(gmm, vs30, measured, z1p0, z2p5):
    gi = GMM_INPUTS
    return plha.get_ground_motion_data(
        gmm,
        vs30,
        measured,
        z1p0,
        z2p5,
        gi["fault_type"],
        gi["rjb"],
        gi["rrup"],
        gi["rx"],
        gi["rx1"],
        gi["ry0"],
        gi["m"],
        gi["ztor"],
        gi["zbor"],
        gi["dip"],
    )


GMM_CASES = [
    ("bssa14", 300.0, False, None, None),
    ("bssa14", 1600.0, False, None, None),
    ("cb14", 300.0, False, None, None),
    ("cb14", 1000.0, False, None, 0.5),
    ("cb14", 300.0, False, None, 4.0),
    ("cb14", 300.0, False, None, 2.0),
    ("cy14", 300.0, True, 0.5, None),
    ("cy14", 1000.0, False, None, None),
    ("ask14", 180.0, True, 0.5, None),
    ("ask14", 250.0, False, None, None),
    ("ask14", 400.0, False, None, None),
    ("ask14", 2000.0, False, None, None),
]


@pytest.mark.parametrize("gmm,vs30,measured,z1p0,z2p5", GMM_CASES)
def test_get_ground_motion_data(gmm, vs30, measured, z1p0, z2p5):
    mu, sigma = _call(gmm, vs30, measured, z1p0, z2p5)
    n = len(GMM_INPUTS["m"])
    assert mu.shape == (n,)
    assert sigma.shape == (n,)
    assert np.all(np.isfinite(mu))
    assert np.all(sigma > 0.0)


def test_ask14_without_z1p0_kwarg():
    gi = GMM_INPUTS
    mu, sigma = ask14.get_im(
        400.0,
        gi["rrup"],
        gi["rx"],
        gi["rx1"],
        gi["ry0"],
        gi["m"],
        gi["fault_type"],
        False,
        gi["dip"],
        gi["ztor"],
    )
    assert np.all(np.isfinite(mu))
    assert np.all(sigma > 0.0)


def test_get_ground_motion_data_invalid_model():
    with pytest.raises(UnboundLocalError):
        _call("not_a_model", 300.0, False, None, None)
