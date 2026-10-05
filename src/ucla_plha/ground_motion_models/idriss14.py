import warnings

import numpy as np


def get_im(vs30, rrup, m, fault_type):
    """
    Output an array of mean and standard deviation values of the natural log of peak horizontal acceleration
    from Idriss (2014), "An NGA-West2 empirical model for estimating the horizontal spectral values generated
    by shallow crustal earthquakes", Earthquake Spectra 30(3), 1155-1177.
    vs30 = time-averaged shear wave velocity in upper 30m [m/s]
    rrup = rupture distance [km]
    m = moment magnitude
    fault_type = style of faulting based on rake. 1 = reverse, 2 = normal, 3 = strike slip

    Notes:
    The model is recommended for vs30 from 450 to 1200 m/s, rrup up to 150 km, and m of 5 or larger. A warning is
    issued for vs30 outside the recommended range, but the model is still evaluated, as in pygmm.
    Normal faulting earthquakes are treated as strike slip, so the style of faulting term applies only to
    reverse earthquakes.
    Coefficients for PGA are from pygmm (idriss_2014-small.csv and idriss_2014-large.csv), which uses
    the coefficients for T = 0.01 s for PGA.
    """
    if (vs30 < 450.0) or (vs30 > 1200.0):
        warnings.warn(
            f"vs30 = {vs30} m/s is outside the 450 to 1200 m/s range recommended for idriss14",
            UserWarning,
            stacklevel=2,
        )

    # Coefficients for m <= 6.75
    alpha_1_small = 7.0887
    alpha_2_small = 0.2058
    # Coefficients for m > 6.75
    alpha_1_large = 9.0138
    alpha_2_large = -0.0794
    # Coefficients that do not depend on magnitude
    alpha_3 = 0.0589
    beta_1 = 2.9935
    beta_2 = -0.2287
    epsilon = -0.854
    gamma = -0.0027
    phi = 0.08
    period = 0.01

    alpha_1 = np.where(m <= 6.75, alpha_1_small, alpha_1_large)
    alpha_2 = np.where(m <= 6.75, alpha_2_small, alpha_2_large)
    frv = np.zeros(len(m), dtype=float)
    frv[fault_type == 1] = 1.0

    mu = (
        alpha_1
        + alpha_2 * m
        + alpha_3 * (8.5 - m) ** 2
        - (beta_1 + beta_2 * m) * np.log(rrup + 10.0)
        + gamma * rrup
        + epsilon * np.log(vs30)
        + phi * frv
    )
    sigma = (
        1.18 + 0.035 * np.log(np.clip(period, 0.05, 3.0)) - 0.06 * np.clip(m, 5.0, 7.5)
    )

    return (mu, sigma)
