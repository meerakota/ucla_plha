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
    delta_n160 = np.exp(1.63 + 9.7/(fc + 0.01) - (15.7 / (fc + 0.01))**2)
    n160cs = n160 + delta_n160
    c_sigma = 1.0 / (18.9 - 2.55 * np.sqrt(n160cs))
    c_sigma = np.min([0.3, c_sigma])
    k_sigma = 1.0 - c_sigma * np.log(sigmavp / pa)
    msf = 6.9 * np.exp(-m / 4.0) - 0.058
    k_sigma = np.min([1.1, k_sigma])
    msf[msf>1.8] = 1.8
    n160cs = np.full(len(m), n160cs)
    k_sigma = np.full(len(m), k_sigma)
    mu_ln_crr = n160cs / 14.1 + (n160cs / 126)**2.0 - (n160cs / 23.6)**3.0 + (n160cs / 25.4)**4 - 2.67 + np.log(k_sigma) + np.log(msf)
    sigma_ln_crr = np.full(len(m), 0.13)
    
    return mu_ln_crr, sigma_ln_crr

def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, d):
    '''
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress, scalar
    sigmavp = vertical effective stress, scalar
    d = depth [m]

    Outputs:
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M
    sigma_ln_csr = standard deviation of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M

    Notes:
    N = number of earthquake events
    '''
    alpha = -1.012 - 1.126 * np.sin(d/11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(d/11.28 + 5.142)
    rd = np.exp(alpha + beta * m)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga
    
    return mu_ln_csr, sigma_ln_csr

def get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, d, n160, fc, fsl, pa):
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
    mu_ln_csr, sigma_ln_csr = get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, d)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    std_ln_fsl = np.sqrt(sigma_ln_csr**2 + sigma_ln_crr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / std_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
