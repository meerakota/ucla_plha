import numpy as np
from scipy.special import ndtr


def get_ln_crr(m, qc1ncs, sigmavp, pa):
    """
    Inputs:
    m = magnitude, Numpy array, dtype=float, length = N
    qc1ncs = scalar clean sand tip resistance
    sigmavp = vertical effective stress, scalar
    pa = atmospheric pressure in same units as sigmavp, scalar

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N

    Notes:
    N = number of earthquake events.
    """
    c_sigma = 1.0 / (37.3 - 8.27 * (qc1ncs) ** 0.264)
    c_sigma = np.minimum(c_sigma, 0.3)
    k_sigma = 1.0 - c_sigma * np.log(sigmavp / pa)
    k_sigma = np.minimum(k_sigma, 1.1)
    msf_max = 1.09 + (qc1ncs / 180) ** 3
    msf_max = np.minimum(msf_max, 2.2)
    msf = 1 + (msf_max - 1) * (8.64 * np.exp(-m / 4) - 1.325)
    mu_ln_crr = (
        qc1ncs / 113
        + (qc1ncs / 1000) ** 2.0
        - (qc1ncs / 140) ** 3.0
        + (qc1ncs / 137) ** 4
        - 2.60
        + np.log(k_sigma)
        + np.log(msf)
    )
    sigma_ln_crr = np.full(len(m), 0.20)

    return mu_ln_crr, sigma_ln_crr


def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth):
    """
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress, scalar
    sigmavp = vertical effective stress, scalar
    depth = depth [m]

    Outputs:
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M
    sigma_ln_csr = standard deviation of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M

    Notes:
    N = number of earthquake events
    """
    alpha = -1.012 - 1.126 * np.sin(depth / 11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(depth / 11.28 + 5.142)
    rd = np.exp(alpha + beta * m)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr


def get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, qc1ncs, fsl, pa):
    """
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress at center of layer, scalar
    sigmavp = vertical effective stress, Numpy array, scalar
    depth = depth [m]
    n160 = energy- and overburden-corrected standard penetration test blow count [-], scalar
    fc = fines content [%], scalar
    fsl = factor of safety values for which we want to compute liquefaction hazard, Numpy array, dtype=float, length = L
    pa = atmospheric pressure in same units as sigmav and sigmavp, scalar

    Outputs:
    fsl_cdfs = cumulative distribution functions for factor of safety against profile manfiestation, Numpy ndarray, dtype=float, shape = N x L
    eps = epsilon for profile manifestation, Numpy ndarray, dtype=float, shape = N x L
    """
    mu_ln_crr, sigma_ln_crr = get_ln_crr(m, qc1ncs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    std_ln_fsl = np.sqrt(sigma_ln_csr**2 + sigma_ln_crr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / std_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
