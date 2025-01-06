import numpy as np
import scipy as sp
import ngl_tools.smt as smt

def pdf(x, mu, sigma):
    return 1.0 / np.sqrt(2.0*np.pi) / sigma * np.exp(-0.5 * ((x - mu)/sigma)**2)

def process_cpt(depth, qt, fs, dGWT, **kwargs):
    gammaw = kwargs.get('gammaw', 9.81)
    Ksat = np.zeros(len(depth))
    Ksat[depth > dGWT] = 1.0
    if('gamma' in kwargs):
        gamma = kwargs.get('gamma')
        sigmav = gamma * depth
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
    

def get_fsl_hazards(ztop, zbot, qc1Ncs_lay, Ic_lay, sigmav_lay, sigmavp_lay, Ksat_lay, mu_ln_pga, sigma_ln_pga, m, rate, fsl, N=50):
    # Define model constants
    lambda_csr = -0.365
    tc = 2.0

    # compute probabilities
    # Prepare variables that can be computed outside of the loop
    crr_hat_lay = smt.get_crr_hat(qc1Ncs_lay)
    log_fsl = np.log(fsl)
    log_min = np.min(log_fsl) - 0.5 * (np.max(log_fsl) - np.min(log_fsl))
    log_max = np.max(log_fsl) + 0.5 * (np.max(log_fsl) - np.min(log_fsl))
    amax = np.logspace(log_min / np.log(10.0), log_max / np.log(10.0), N)
    log_amax = np.log(amax)
    dx = log_amax[1] - log_amax[0]
    log_amax_diff = log_amax[:, np.newaxis] - log_amax
    pfs = smt.get_pfs(Ic_lay)
    pfmt = smt.get_pfmt(ztop, Ic_lay)
    t = zbot - ztop

    # Initialize result array and perform convolution numerically via matrix multiplication
    result = np.zeros((len(m), len(amax)))
    for i in range(len(m)):
        pdf_array = pdf(log_amax_diff, mu_ln_pga[i], sigma_ln_pga[i])
        w = dx * pdf_array
        csrm_lay = smt.get_csrm(amax[:, np.newaxis], m[i], sigmav_lay, sigmavp_lay, 0.5 * (ztop + zbot), qc1Ncs_lay)
        csrm_hat_lay = smt.box_cox(csrm_lay, lambda_csr)
        pfts = smt.get_pfts(csrm_hat_lay, crr_hat_lay, Ksat_lay)
        pmp = 1.0 - np.prod((1.0 - pfmt * pfts * pfs) ** (t / tc), axis=1)
        temp = pmp @ w
        result[i] = sp.interpolate.pchip_interpolate(amax, temp.T, fsl, axis=0)
    fsl_hazards = rate[:, np.newaxis] * result
    return fsl_hazards