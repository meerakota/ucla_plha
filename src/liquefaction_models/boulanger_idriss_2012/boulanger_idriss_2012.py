#Boulanger and Idriss 2012
#https://doi.org/10.1061/(ASCE)GT.1943-5606.0000700
import numpy as np
from scipy.special import ndtr

def get_crr(m, n160, fc, sigmavp):
    pa = 101.325
    delta_n160 = np.exp(1.63 + 9.7/(fc + 0.01) - (15.7 / (fc + 0.01))**2)
    n160cs = n160 + delta_n160
    c_sigma = 1.0 / (18.9 - 2.55 * np.sqrt(n160cs))
    k_sigma = 1.0 - c_sigma * np.log(sigmavp / pa)
    msf = 6.9 * np.exp(-m / 4.0) - 0.058
    c_sigma = np.min([0.3, c_sigma])
    k_sigma = np.min([1.1, k_sigma])
    msf[msf>1.8] = 1.8
    n160cs = np.full(len(m), n160cs)
    k_sigma = np.full(len(m), k_sigma)
    mu_crr = n160cs / 14.1 + (n160cs / 126)**2.0 - (n160cs / 23.6)**3.0 + (n160cs / 25.4)**4 - 2.67 + np.log(k_sigma) + np.log(msf)
    std_crr = np.full(len(m), 0.13)
    return [mu_crr, std_crr]

def get_csr(mu_pga, std_pga, m, sigmav, sigmavp, d):
    alpha = -1.012 - 1.126 * np.sin(d/11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(d/11.28 + 5.142)
    rd = np.exp(alpha + beta * m)
    mu_csr = mu_pga + np.log(0.65 * sigmav / sigmavp * rd)
    std_csr = std_pga
    return [mu_csr, std_csr]

def get_fsl_hazards(mu_pga, std_pga, rate, fsl, m, sigmav, sigmavp, d, n160, fc, pa):
    delta_n160 = np.exp(1.63 + 9.7/(fc + 0.01) - (15.7 / (fc + 0.01))**2)
    n160cs = n160 + delta_n160
    c_sigma = 1.0 / (18.9 - 2.55 * np.sqrt(n160cs))
    k_sigma = 1.0 - c_sigma * np.log(sigmavp / pa)
    msf = 6.9 * np.exp(-m / 4.0) - 0.058
    c_sigma = np.min([0.3, c_sigma])
    k_sigma = np.min([1.1, k_sigma])
    msf[msf>1.8] = 1.8
    n160cs = np.full(len(m), n160cs)
    k_sigma = np.full(len(m), k_sigma)
    mu_crr = n160cs / 14.1 + (n160cs / 126)**2.0 - (n160cs / 23.6)**3.0 + (n160cs / 25.4)**4 - 2.67 + np.log(k_sigma) + np.log(msf)
    std_crr = np.full(len(m), 0.13)
    alpha = -1.012 - 1.126 * np.sin(d/11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(d/11.28 + 5.142)
    rd = np.exp(alpha + beta * m)
    mu_csr = mu_pga + np.log(0.65 * sigmav / sigmavp * rd)
    std_csr = std_pga
    mu_fsl = np.log(mu_csr) - np.log(mu_crr)
    std_fsl = np.sqrt(std_csr**2 + std_crr**2)
    fsl_hazards = ndtr((np.log(fsl) - mu_fsl[:, np.newaxis]) / std_fsl[:, np.newaxis]) * rate[:, np.newaxis]
    return fsl_hazards
