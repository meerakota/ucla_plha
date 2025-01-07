import numpy as np
from scipy.special import ndtr

def get_crr(m, n160, fc, sigmavp):
    pa = 101.325
    mu_crr = (n160*(1+0.00167*fc) - 27.352 * np.log(m) - 3.958 * np.log(sigmavp / pa) + 0.089 * fc + 16.084)/11.771
    std_crr = np.full(len(m), 2.95 / 11.771)
    return [mu_crr, std_crr]

def get_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, vs12, d):
    rd_num = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (-d + 0.0785 * vs12 + 7.586)))
    rd_den = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    rd = rd_num / rd_den
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga
    return [mu_ln_csr, sigma_ln_csr]

def get_fsl_hazards(mu_ln_pga, sigma_ln_pga, rate, fsl, m, sigmav, sigmavp, vs12, d, n160, fc, pa):
    mu_ln_crr = (n160*(1+0.00167*fc) - 27.352 * np.log(m) - 3.958 * np.log(sigmavp / pa) + 0.089 * fc + 16.084)/11.771
    sigma_ln_crr = np.full(len(m), 2.95 / 11.771)
    rd_num = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (-d + 0.0785 * vs12 + 7.586)))
    rd_den = 1.0 + (-23.013 - 2.949 * np.exp(mu_ln_pga) + 0.999 * m + 0.0525 * vs12) / (16.258 + 0.201 * np.exp(0.341 * (0.0785 * vs12 + 7.586)))
    rd = rd_num / rd_den
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga
    mu_ln_fsl = mu_ln_csr - mu_ln_crr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    fsl_hazards = ndtr((np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]) * rate[:, np.newaxis]
    return fsl_hazards