import numpy as np
from scipy.special import ndtr

def get_ln_crr(m, n160, fc, sigmavp, pa):
    '''
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
    '''
    mu_ln_crr = (n160*(1+0.00167*fc) - 27.352 * np.log(m) - 3.958 * np.log(sigmavp / pa) + 0.089 * fc + 16.084)/11.771
    sigma_ln_crr = np.full(len(m), 2.95 / 11.771)
    
    return mu_ln_crr, sigma_ln_crr

def get_rd(mu_ln_pga, m, vs12, d):
    '''
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
    '''
    if(d < 20):
        rd_num = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (-d + 0.0785 * vs12 + 7.586)))
        rd_den = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    else:
        rd_num = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (-20 + 0.0785 * vs12 + 7.586)))
        rd_den = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    rd = rd_num / rd_den
    
    return rd

def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d):
    '''
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
    '''
    rd = get_rd(mu_ln_pga, m, vs12, d)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr

def get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d, n160, fc, fsl, pa):
    '''
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

    Outputs:
    fsl_cdfs = cumulative distribution functions for factor of safety against profile manfiestation, Numpy ndarray, dtype=float, shape = N x L
    eps = epsilon for profile manifestation, Numpy ndarray, dtype=float, shape = N x L
    '''
    mu_ln_crr, sigma_ln_crr = get_ln_crr(m, n160, fc, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
