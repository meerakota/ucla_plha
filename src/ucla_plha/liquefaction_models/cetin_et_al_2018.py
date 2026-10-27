"""Cetin et al. (2018) SPT-based liquefaction triggering model.

References:
    Cetin, K.O., Seed, R.B., Kayen, R.E., Moss, R.E.S., Bilge, H.T., Ilgac, M., and Chowdhury, K.
    (2018a). SPT-based probabilistic and deterministic assessment of seismic soil liquefaction
    triggering hazard. Soil Dynamics and Earthquake Engineering, 115, 698-709.

    Cetin, K.O., Seed, R.B., Kayen, R.E., Moss, R.E.S., Bilge, H.T., Ilgac, M., and Chowdhury, K.
    (2018b). The use of the SPT-based seismic soil liquefaction triggering evaluation methodology
    in engineering hazard assessments. MethodsX, 5, 1556-1575. doi:10.1016/j.mex.2018.11.016
"""
import numpy as np
from scipy.special import ndtr


# Coefficients of Cetin et al. (2018a), Table 1 of Cetin et al. (2018b)
THETA_1 = 0.00167
THETA_2 = 27.352
THETA_3 = 3.958
THETA_4 = 0.089
THETA_5 = 16.084
THETA_6 = 11.771
SIGMA_EPSILON = 2.95
# Scaling factor of the input parameter uncertainty, Eq. 4 of Cetin et al. (2018b)
THETA_7 = 0.392


def get_sigma_total(n160, fc, parameter_uncertainty=False, sigma_n160=0.0, sigma_fc=0.0,
                    cov_sigmavp=0.0, cov_csr=0.0):
    """
    Standard deviation of the limit state function (in the units of N1,60), sigma_tot.

    Without parameter uncertainty, sigma_tot is the model uncertainty sigma_epsilon = 2.95.
    With parameter uncertainty, sigma_tot**2 = theta_7 * sigma_input**2 + sigma_epsilon**2
    (Eq. 4 of Cetin et al. 2018b), where (Eq. 5)

        sigma_input**2 = theta_6**2 * cov_csr**2 + sigma_n160**2 * (1 + theta_1 * fc)**2
                         + sigma_fc**2 * (theta_1 * n160 + theta_4)**2
                         + theta_3**2 * cov_sigmavp**2

    Inputs:
    n160 = energy- and overburden-corrected standard penetration test blow count [-], scalar
    fc = fines content [%], scalar
    parameter_uncertainty = include the input parameter uncertainty, bool
    sigma_n160 = standard deviation of n160 [-], scalar
    sigma_fc = standard deviation of fc [%], scalar
    cov_sigmavp = coefficient of variation of the vertical effective stress [-], scalar
    cov_csr = coefficient of variation of the cyclic stress ratio [-], scalar. In ucla_plha the
        variability of the ground motion is integrated in the hazard calculation, so cov_csr
        should only represent the uncertainty in estimating CSR for a given ground motion
        (e.g., rd and the stresses), not the ground motion variability.

    Outputs:
    sigma_tot = standard deviation of the limit state function, scalar
    """
    if not parameter_uncertainty:
        return SIGMA_EPSILON
    sigma_input_sq = (
        THETA_6**2 * cov_csr**2
        + sigma_n160**2 * (1 + THETA_1 * fc) ** 2
        + sigma_fc**2 * (THETA_1 * n160 + THETA_4) ** 2
        + THETA_3**2 * cov_sigmavp**2
    )
    return np.sqrt(THETA_7 * sigma_input_sq + SIGMA_EPSILON**2)


def get_ln_crr(m, n160, fc, sigmavp, pa, parameter_uncertainty=False, sigma_n160=0.0,
               sigma_fc=0.0, cov_sigmavp=0.0, cov_csr=0.0):
    """
    Inputs:
    m = magnitude, Numpy array, dtype=float, length = N
    n160 = energy- and overburden-corrected standard penetration test blow count [-], scalar
    fc = fines content [%], scalar
    sigmavp = vertical effective stress, scalar
    pa = atmospheric pressure in same units as sigmavp, scalar
    parameter_uncertainty, sigma_n160, sigma_fc, cov_sigmavp, cov_csr = input parameter
        uncertainty (see get_sigma_total); by default only the model uncertainty is included

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N

    Notes:
    N = number of earthquake events.
    """
    mu_ln_crr = (
        n160 * (1 + THETA_1 * fc)
        - THETA_2 * np.log(m)
        - THETA_3 * np.log(sigmavp / pa)
        + THETA_4 * fc
        + THETA_5
    ) / THETA_6
    sigma_tot = get_sigma_total(
        n160, fc, parameter_uncertainty, sigma_n160, sigma_fc, cov_sigmavp, cov_csr
    )
    sigma_ln_crr = np.full(len(m), sigma_tot / THETA_6)

    return mu_ln_crr, sigma_ln_crr


def get_rd(mu_ln_pga, m, vs12, d):
    """
    Inputs:
    m = magnitude, Numpy array, dtype=float, length = N
    n160 = energy- and overburden-corrected standard penetration test blow count [-], scalar
    fc = fines content [%], scalar
    sigmavp = vertical effective stress, scalar
    pa = atmospheric pressure in same units as sigmavp, scalar

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N

    Notes:
    N = number of earthquake events.
    """
    if d < 20:
        rd_num = 1.0 + (
            -23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12
        ) / (16.258 + 0.201 * np.exp(0.341 * (-d + 0.0785 * vs12 + 7.586)))
        rd_den = 1.0 + (
            -23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12
        ) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    else:
        rd_num = 1.0 + (
            -23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12
        ) / (16.258 + 0.201 * np.exp(0.341 * (-20 + 0.0785 * vs12 + 7.586)))
        rd_den = 1.0 + (
            -23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12
        ) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    rd = rd_num / rd_den

    return rd


def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d):
    """
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress, scalar
    sigmavp = vertical effective stress, scalar
    vs12 = time-averaged shear wave velocity in upper 12 m [m/s], scalar
    d = depth [m]

    Outputs:
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M
    sigma_ln_csr = standard deviation of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M

    Notes:
    N = number of earthquake events
    """
    rd = get_rd(mu_ln_pga, m, vs12, d)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr


def get_fsl_cdfs(
    mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d, n160, fc, fsl, pa,
    parameter_uncertainty=False, sigma_n160=0.0, sigma_fc=0.0, cov_sigmavp=0.0, cov_csr=0.0,
):
    """
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress at center of layer, scalar
    sigmavp = vertical effective stress, Numpy array, scalar
    d = depth [m]
    n160 = energy- and overburden-corrected standard penetration test blow count [-], scalar
    fc = fines content [%], scalar
    fsl = factor of safety values for which we want to compute liquefaction hazard, Numpy array, dtype=float, length = L
    pa = atmospheric pressure in same units as sigmav and sigmavp, scalar
    parameter_uncertainty, sigma_n160, sigma_fc, cov_sigmavp, cov_csr = input parameter
        uncertainty of Cetin et al. (2018b), see get_sigma_total. By default (False), only the
        model uncertainty is included, as before.

    Outputs:
    fsl_cdfs = cumulative distribution functions for factor of safety against profile manfiestation, Numpy ndarray, dtype=float, shape = N x L
    eps = epsilon for profile manifestation, Numpy ndarray, dtype=float, shape = N x L
    """
    mu_ln_crr, sigma_ln_crr = get_ln_crr(
        m, n160, fc, sigmavp, pa, parameter_uncertainty, sigma_n160, sigma_fc, cov_sigmavp,
        cov_csr,
    )
    mu_ln_csr, sigma_ln_csr = get_ln_csr(
        mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d
    )
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
