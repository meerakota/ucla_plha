import numpy as np
import scipy as sp
import ngl_tools.smt as smt
from scipy.stats import norm
from scipy.special import ndtr

def pdf(x, mu, sigma):
    return 1.0 / np.sqrt(2.0*np.pi) / sigma * np.exp(-0.5 * ((x - mu)/sigma)**2)

def process_cpt(depth, qt, fs, dGWT, **kwargs):
    gammaw = kwargs.get('gammaw', 9.81)
    Ksat = np.zeros(len(depth))
    Ksat[depth > dGWT] = 1.0
    if('gamma' in kwargs):
        gamma = kwargs.get('gamma')
        sigmav = gamma * depth
        if(dGWT < 0):
            sigmav -= dGWT * gammaw
        u = gammaw * (depth - dGWT)
        sigmavp = sigmav - u
    elif('sigmav' in kwargs):
        sigmav = kwargs.get('sigmav')
        sigmavp = kwargs.get('sigmavp')
    # apply inverse-filter algorithm
    qt_inv, fs_inv, Ic_inv = smt.cpt_inverse_filter(qt, depth, fs=fs, sigmav=sigmav, sigmavp=sigmavp, low_pass=True, smooth=True, remove_interface=True)
    Ic_inv, Qtn_inv, Fr_inv = smt.get_Ic_Qtn_Fr(qt_inv, fs_inv, sigmav, sigmavp)

    # Compute fines content and overburden- and fines-corrected tip resistance for inverse-filtered CPT data
    FC = smt.get_FC_from_Ic(Ic_inv, 0.0)
    qc1N_inv, qc1Ncs_inv = smt.get_qc1N_qc1Ncs(qt_inv, fs_inv, sigmav, sigmavp, FC)

    # apply layering algorithm
    ztop, zbot, qc1Ncs_lay, Ic_lay = smt.cpt_layering(qc1Ncs_inv, Ic_inv, depth, dGWT=dGWT, Nmin=1, Nmax=None)
    # insert layer break at groundwater table depth if needed
    if((dGWT not in ztop) and (dGWT not in zbot) and (dGWT > np.min(ztop)) and (dGWT < np.max(zbot))):
        Ic_dGWT = Ic_lay[ztop<dGWT][-1]
        qc1Ncs_dGWT = qc1Ncs_lay[ztop<dGWT][-1]
        Ic_lay = np.hstack((Ic_lay[zbot<dGWT], Ic_dGWT, Ic_lay[zbot>dGWT]))
        qc1Ncs_lay = np.hstack((qc1Ncs_lay[zbot<dGWT], qc1Ncs_dGWT, qc1Ncs_lay[zbot>dGWT]))
        ztop = np.hstack((ztop[ztop<dGWT], dGWT, ztop[ztop>dGWT]))
        zbot = np.hstack((zbot[zbot<dGWT], dGWT, zbot[zbot>dGWT]))
    sigmav_lay = np.interp(0.5 * (ztop + zbot), depth, sigmav)
    sigmavp_lay = np.interp(0.5 * (ztop + zbot), depth, sigmavp)
    Ksat_lay = np.interp(0.5 * (ztop + zbot), depth, Ksat)
    return (ztop, zbot, qc1Ncs_lay, Ic_lay, sigmav_lay, sigmavp_lay, Ksat_lay)

def get_ln_crr(mu_ln_pga, m, qc1Ncs, sigmavp, pa):
    '''
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    qc1Ncs = overburden- and fines-corrected cone tip resistance [-], Numpy array, dtype=float, length = M
    sigmavp = vertical effective stress, Numpy array, dtype=float, length = M
    pa = atmospheric pressure in same units as sigmavp, scalar

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N x M
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N x M

    Notes:
    N = number of earthquake events.
    M = number of soil layers in profile.
    We are using the mean value of pga to compute msf, but pga is actually a random variable.
    '''
    lambda_dr = 1.2259
    dr = 47.8 * qc1Ncs ** 0.264 - 106.3
    dr[dr < 0] = 0
    dr[dr > 160] = 160
    dr_hat = (dr ** lambda_dr - 1) / lambda_dr
    ksigma = (sigmavp / pa) ** (-0.0015 * dr)
    ksigma[ksigma > 1.2] = 1.2
    neq = np.exp(0.46305 - 0.4082 * np.exp(mu_ln_pga) + 0.2332 * m)
    msf = (14 / neq) ** 0.2
    mu_ln_crr = -3.204 + 0.011 * dr_hat[np.newaxis, :] + np.log(ksigma[np.newaxis, :]) + np.log(msf[:, np.newaxis])
    sigma_ln_crr = 0.342

    return mu_ln_crr, sigma_ln_crr

def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, sigmav, sigmavp):
    '''
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    ztop = depth to top of layer, Numpy array, dtype=float, length = M
    zbot = depth to bottom of layer, Numpy array, dtype=float, length = M
    sigmav = vertical total stress at center of layer, Numpy array, dtype=float, length = M
    sigmavp = vertical effective stress at center of layer, Numpy array, dtype=float, length = M
    pa = atmospheric pressure in same units as sigmav and sigmavp, scalar

    Outputs:
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M
    sigma_ln_csr = standard deviation of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M

    Notes:
    N = number of earthquake events
    M = number of soil layers in profile
    We have include magnitude scale factor here because it depends on magnitude and peak acceleration. But we have moved ksigma
    to get_ln_crr because it depends on relative density. This avoids repetition.
    '''
    d = (ztop + zbot) / 2.0
    alpha = np.exp(-4.373 + 0.4491 * m)
    beta = -20.11 + 6.247 * m
    rd = (1.0 - alpha[:, np.newaxis]) * np.exp(-d / beta[:, np.newaxis]) + alpha[:, np.newaxis]
    mu_ln_csr = mu_ln_pga[:, np.newaxis] + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr

def get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa):
    '''
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    ztop = depth to top of layer, Numpy array, dtype=float, length = M
    zbot = depth to bottom of layer, Numpy array, dtype=float, length = M
    qc1Ncs = overburden- and fines-corrected cone tip resistance [-], Numpy array, dtype=float, length = M
    Ic = soil behavior type index [-], Numpy array, dtype=float, length = M
    sigmav = vertical total stress at center of layer, Numpy array, dtype=float, length = M
    sigmavp = vertical effective stress, Numpy array, dtype=float, length = M
    Ksat = soil saturation term (1 below groundwater, 0 above), Numpy array, dtype=float, length = M
    fsl = factor of safety values for which we want to compute liquefaction hazard, Numpy array, dtype=float, length = L
    pa = atmospheric pressure in same units as sigmav and sigmavp, scalar

    Outputs:
    fsl_cdfs = cumulative distribution functions for factor of safety against profile manfiestation, Numpy ndarray, dtype=float, shape = N x L
    eps = epsilon for profile manifestation, Numpy ndarray, dtype=float, shape = N x L

    Notes:
    N = number of earthquake events
    M = number of soil layers in profile
    L = number of elements in factor of safety array
    epsilon is generally defined as the number of standard deviations from the mean, and is generally applied to normally
    distributed random variables. Commonly, a unit normal cumulative distribution operator is used to compute the cumulative probability
    density associated with a particular epsilon value (e.g., +1 epsilon = 0.84). In this case, the probability of profile manifestation 
    is not normally distributed, and the functional form of the distribution is not directly known. We are sampling the CDF of 
    manifestation at specified factor of safety values, but we are not evaluating the functional form of the distribution. Therefore, 
    cumulative probability density cannot be computed directly from the mean and standard deviation. We made a modeling decision to 
    compute epsilon as the inverse unit normal cumulative distribution operator of the computed CDF values. In this manner, users 
    can interpret the reported epsilon values as corresponding to specific probabilities of exceedance (e.g., +1 epsilon = 0.84), 
    but should not interpret them as being equal to the number of standard deviations from the mean. Furthermore, the reported epsilons
    correspond to manifestation rather than triggering, and therefore include susceptibility, manifestation conditioned on triggering, 
    and saturation terms. Cumulative distributions of manifestation do not necessarily reach 1.0 at high factor of safety because,
    for example, a soil layer may have a low probability of susceptibility or manifestation conditioned on triggering. Manifestation is 
    therefore a mixed random variable with some probability mass lumped at fsl = infinity. This tends to bias the disaggregation toward 
    low epsilons because we don't plot it to fsl = infinity.
    '''
    tc = 2.0
    mu_ln_crr, sigma_ln_crr = get_ln_crr(mu_ln_pga, m, qc1Ncs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = get_ln_csr(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, sigmav, sigmavp)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr ** 2 + sigma_ln_csr ** 2)
    eps_triggering = (np.log(fsl[np.newaxis, np.newaxis, :]) - mu_ln_fsl[:,:,np.newaxis]) / sigma_ln_fsl[:, np.newaxis, np.newaxis]
    pfts = ndtr(eps_triggering)
    pfs = 1.0 - 1.0 / (1.0 + np.exp(-14.8 * (Ic / 2.635 - 1.0)))  
    pfmt = 1.0 / (1.0 + np.exp(-(7.613 - 0.338 * ztop - 3.042 * Ic))) 
    t = zbot - ztop
    fsl_cdfs = 1.0 - np.prod((1.0 - pfmt[np.newaxis, :, np.newaxis] * pfts * pfs[np.newaxis, :, np.newaxis] * Ksat[np.newaxis, :, np.newaxis]) ** (t[np.newaxis, :, np.newaxis] / tc), axis=1)
    eps_manifestation = norm.ppf(fsl_cdfs)

    return fsl_cdfs, eps_manifestation
