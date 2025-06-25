import numpy as np
from scipy.special import ndtr

def get_ln_crr(m, qc, fs, sigmavp, pa):
    '''
    Inputs:
    m = magnitude, Numpy array, dtype=float, length = N
    qc = cone penetration test tip resistance in same units as pa, scalar
    fs = cone penetration test sleeve friction in same units as pa, scalar
    sigmavp = vertical effective stress in same units as pa, scalar
    pa = atmospheric pressure, scalar

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N

    Notes:
    N = number of earthquake events.
    '''
    rf = fs / qc
    f1 = 0.78 * qc ** -0.33
    f2 = -(-0.32 * qc ** -0.35 + 0.49)
    f3 = np.abs(np.log10(10.0 + qc)) ** 1.21
    c = f1 * (rf / f3) ** f2
    Cq = (pa / sigmavp) ** c
    qc1 = Cq * qc
    mu_ln_crr = qc1 ** 1.045 + qc1 * (0.110 * rf) + (0.001 * rf) + c * (1.0 + 0.850 * rf) - 0.848 * np.log(m) - 0.002 * np.log(sigmavp) - 20.923
    sigma_ln_crr = np.full(len(m), 1.632 / 7.177)
    
    return mu_ln_crr, sigma_ln_crr

def get_rd(mu_ln_pga, m, d):
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
    amax = np.exp(mu_ln_pga)
    if(d < 20):
        rd_num = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (10.567 + 0.089 * np.exp(0.089 * (-d * 3.28 - 7.760 * amax + 78.576)))
        rd_den = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (10.567 + 0.089 * np.exp(0.089 * (- 7.760 * amax + 78.576)))
        rd = rd_num / rd_den
    else:
        rd_num = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (10.567 + 0.089 * np.exp(0.089 * (-d * 3.28 - 7.760 * amax + 78.576)))
        rd_den = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (10.567 + 0.089 * np.exp(0.089 * (- 7.760 * amax + 78.576)))
        rd = rd_num / rd_den - 0.0014 * (d * 3.28 - 65.0)
    
    return rd

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
    m_adjusted = np.copy(m)
    m_adjusted[m_adjusted < 5.5] = 5.5
    m_adjusted[m_adjusted > 8.5] = 8.5
    dwfm = 17.84 * m_adjusted ** -1.43    
    rd = get_rd(mu_ln_pga, m, d)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd) - np.log(dwfm)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr

def get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, d, qc1, fs, fsl, pa):
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
    mu_ln_crr, sigma_ln_crr = get_ln_crr(m, qc1, fs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, d)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
