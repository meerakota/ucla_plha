import numpy as np
from scipy.special import ndtr, ndtri

from ucla_plha import pga_grid

# Default upper limit of amax [g] in rd for the probability conditioned on PGA
# (get_conditional_fsl_cdfs, config "rd_amax_cap"). The rd of Cetin et al. (2004) used by Moss et
# al. (2006) decreases with amax to the lower limit 0.2 at about 2.5 to 3.5 g and jumps to the upper
# limit 1.0 at about 3.5 to 4 g, where its denominator changes sign; with the cap, rd above 1 g is
# rd at 1 g.
RD_AMAX_CAP = 1.0

# CSR without duration weighting: Moss et al. (2006) compare Eq. 20 and 21, which have their own
# magnitude term (-0.848 ln Mw), with the CSR of the event: the note to Eq. 21 says it "produces
# results for the input magnitude CRR_Mw=?, whereas Figs. 6 and 10 give results for a mean magnitude
# of 7.5 (CSR*)", and the CSR of the case histories in their Table 1 is 0.65 amax (sigmav / sigmav')
# rd, not CSR / DWF (e.g. 1987 Edgecumbe, Mw 6.6, and 1981 Westmorland, Mw 5.9). ucla_plha 2.x
# divided the CSR by DWF = 17.84 Mw^-1.43 as well, which double counted the magnitude effect (FS
# higher by the DWF, about 1.4 at Mw 6, 1.1 at Mw 7, and 0.91 at Mw 8). The 2.x Moss et al. (2006)
# results, which divided CSR by the DWF, can't be reproduced with this version.


def get_ln_crr(m, mu_ln_csr, qc, fs, sigmavp, pa, nonnegative_delta_qc=False):
    """
    Inputs:
    m = magnitude, Numpy array, dtype=float, length = N
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy array, dtype=float, length = N
    qc = cone penetration tip resistance in MPa
    fs = sleeve resistance in MPa
    sigmavp = vertical effective stress in kPa
    pa = atmospheric pressure in kPa, 101.325 kPa
    nonnegative_delta_qc = if True, delta_qc = max(0, x1 ln(CSR) + x2) (get_conditional_fsl_cdfs).
        Moss et al. (2006, Eq. 18) bound the "fines" adjustment to Rf of 0.5 to 5% (delta_qc = 0
        at Rf <= 0.5%) and describe it as a suppression of liquefiability; x1 ln(CSR) + x2 < 0 only
        for CSR < exp(-x2 / x1) = exp(-3.84) = 0.021 (for any Rf), outside their data, where it
        can also make qc1 + delta_qc negative (qc1_mod**1.045 is then nan). The CSR conditioned on
        PGA reaches such values in the lower tail of the PGA distribution. False (default): the
        ucla_plha 2.x formula.

    Outputs:
    mu_ln_crr = mean of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N
    sigma_ln_crr = standard deviation of natural logs of cyclic resistance ratio, Numpy array, dtype=float, length = N

    Notes:
    N = number of earthquake events.
    #"""
    rf = fs / qc * 100
    rf = np.min([rf, 5.0])
    rf = np.max([rf, 0.5])
    f1 = 0.78 * qc**-0.33
    f2 = -(-0.32 * qc**-0.35 + 0.49)
    f3 = np.abs(np.log10(10.0 + qc)) ** 1.21
    c = f1 * (rf / f3) ** f2
    Cq = (pa / sigmavp) ** c
    Cq = np.minimum(Cq, 1.7)
    qc1 = Cq * qc
    x1 = 0.38 * np.min([rf, 5.0]) - 0.19
    x2 = 1.46 * np.min([rf, 5.0]) - 0.73
    if rf > 0.5:
        delta_qc = x1 * mu_ln_csr + x2
        if nonnegative_delta_qc:
            delta_qc = np.maximum(delta_qc, 0.0)
    else:
        delta_qc = 0.0
    qc1_mod = qc1 + delta_qc
    mu_ln_crr = (
        qc1_mod**1.045
        + qc1_mod * (0.110 * rf)
        + (0.001 * rf)
        + c * (1.0 + 0.850 * rf)
        - 0.848 * np.log(m)
        - 0.002 * np.log(sigmavp)
        - 20.923
    ) / 7.177
    sigma_ln_crr = np.full(len(m), 1.632 / 7.177)

    return mu_ln_crr, sigma_ln_crr


def get_rd(mu_ln_pga, m, depth, amax_cap=None):
    """
    Inputs:
    mu_ln_pga = natural log of peak acceleration amax [g], Numpy array, dtype=float, length = N
    depth = depth [m], scalar
    amax_cap = upper limit of amax [g] in rd; None (default): no limit (ucla_plha 2.x)

    rd of Cetin et al. (2004), Moss et al. (2006) Eqs. 8 (depth < 20 m) and 9 (depth >= 20 m),
    with depth in feet (3.28 depth). Eq. 9 of Moss et al. (2006) has 78.567 in the exponent where
    Eq. 8 has 78.576, likely a typo; 78.576 is used for both.
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
    amax = np.exp(mu_ln_pga)
    if amax_cap is not None:
        amax = np.minimum(amax, amax_cap)
    if depth < 20:
        rd_num = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (
            10.567 + 0.089 * np.exp(0.089 * (-depth * 3.28 - 7.760 * amax + 78.576))
        )
        rd_den = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (
            10.567 + 0.089 * np.exp(0.089 * (-7.760 * amax + 78.576))
        )
        rd = rd_num / rd_den
    else:
        rd_num = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (
            10.567 + 0.089 * np.exp(0.089 * (-depth * 3.28 - 7.760 * amax + 78.576))
        )
        rd_den = 1.0 + (-9.147 - 4.173 * amax + 0.652 * m) / (
            10.567 + 0.089 * np.exp(0.089 * (-7.760 * amax + 78.576))
        )
        rd = rd_num / rd_den - 0.0014 * (depth * 3.28 - 65.0)

    # For large values of amax, rd can become negative or larger than 1.0. Constrain
    # these values to lie within 0.2 to 1.0, which is the approximate range in Fig. 4
    # of Moss et al. (2006)
    rd[rd < 0.2] = 0.2
    rd[rd > 1.0] = 1.0
    return rd


def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, rd_amax_cap=None):
    """
    Inputs:
    mu_ln_pga = mean of natural logs of peak acceleration [g], Numpy array, dtype=float, length = N
    sigma_ln_pga = standard deviation of natural logs of peak acceleration [-], Numpy array, dtype=float, length = N
    m = magnitude, Numpy array, dtype=float, length = N
    sigmav = vertical total stress, scalar
    sigmavp = vertical effective stress, scalar
    depth = depth [m]
    rd_amax_cap = upper limit of amax [g] in rd (get_rd); None (default): no limit

    Outputs:
    mu_ln_csr = mean of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M
    sigma_ln_csr = standard deviation of natural logs of cyclic stress ratio, Numpy ndarray, dtype=float, shape = N x M

    Notes:
    N = number of earthquake events
    CSR = 0.65 amax (sigmav / sigmav') rd, without duration weighting (see the note at the top of
    this module; ucla_plha 2.x divided it by DWF = 17.84 Mw^-1.43)
    """
    rd = get_rd(mu_ln_pga, m, depth, rd_amax_cap)
    mu_ln_csr = mu_ln_pga + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr


def get_conditional_fsl_cdfs(
    m, ln_pga, sigmav, sigmavp, depth, qc, fs, fsl, pa, rd_amax_cap=RD_AMAX_CAP
):
    """
    Cumulative distribution function of the factor of safety conditioned on magnitude and peak
    acceleration.

    Inputs:
    m = magnitude, Numpy array, dtype=float, length = K
    ln_pga = natural log of peak acceleration [g], Numpy array, dtype=float, length = K
    sigmav, sigmavp, depth, qc, fs, fsl, pa: as for get_fsl_cdfs
    rd_amax_cap = upper limit of amax [g] in rd (default RD_AMAX_CAP = 1 g); None: no limit

    Outputs:
    fsl_cdfs = P(FS < fsl | m, PGA) = Phi((ln fsl - mu_ln_crr(m, CSR) + ln CSR(m, PGA)) / (1.632 / 7.177)),
        with rd (amax = PGA, at most rd_amax_cap) and the Delta qc adjustment of the CRR (x1 ln CSR + x2) evaluated at the
        given PGA (delta_qc not below 0, see get_ln_crr), Numpy ndarray, dtype=float, shape = K x L
    """
    m = np.asarray(m, dtype=float)
    ln_pga = np.asarray(ln_pga, dtype=float)
    mu_ln_csr, _ = get_ln_csr(ln_pga, 0.0, m, sigmav, sigmavp, depth, rd_amax_cap)
    mu_ln_crr, sigma_ln_crr = get_ln_crr(
        m, mu_ln_csr, qc, fs, sigmavp, pa, nonnegative_delta_qc=True
    )
    return ndtr(
        (np.log(fsl)[np.newaxis, :] - (mu_ln_crr - mu_ln_csr)[:, np.newaxis])
        / sigma_ln_crr[:, np.newaxis]
    )


def get_fsl_cdfs(
    mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth, qc, fs, fsl, pa,
    pga_uncertainty="conditional", truncation_level=None, n_quadrature=128,
    rd_amax_cap=RD_AMAX_CAP,
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
    pga_uncertainty = treatment of the uncertainty of the peak acceleration, string:
        "conditional" (default): P(FS < fsl | m, PGA) = Phi((ln fsl - mu_ln_crr + ln CSR(m, PGA)) / sigma_ln_crr)
            (get_conditional_fsl_cdfs) is integrated over the distribution of ln(PGA) of the event, truncated at
            truncation_level standard deviations as in the PSHA (Gauss-Legendre quadrature with n_quadrature points;
            plha.get_hazard integrates it on the magnitude x PGA grid of ucla_plha.pga_grid instead).
        "closed_form": the treatment of ucla_plha 2.x, untruncated: Phi((ln fsl - mu_ln_fsl) / sqrt(sigma_ln_crr^2 + sigma_ln_pga^2))
            with rd and the Delta qc adjustment of the CRR evaluated at the median PGA (and the median CSR), but without
            the duration weighting of 2.x: the 2.x Moss et al. (2006) results, which divided CSR by the DWF, can't be
            reproduced with this version (the DWF double counted the magnitude effect; see the note at the top of this
            module, Eq. 21 note and Table 1 of Moss et al. 2006).
    truncation_level = upper truncation of the distribution of ln(PGA) in standard deviations ("conditional"
        only), as in the PSHA (plha.get_exceedance); None (default) is untruncated
    n_quadrature = number of quadrature points ("conditional")
    rd_amax_cap = upper limit of amax [g] in rd ("conditional" only; default RD_AMAX_CAP = 1 g,
        None: no limit). "closed_form" uses rd of the median PGA without a limit, as ucla_plha 2.x.

    Outputs:
    fsl_cdfs = cumulative distribution functions for factor of safety against profile manfiestation, Numpy ndarray, dtype=float, shape = N x L
    eps = epsilon for profile manifestation, Numpy ndarray, dtype=float, shape = N x L ("conditional": the
        inverse standard normal CDF of fsl_cdfs)
    """
    if pga_uncertainty == "conditional":

        def conditional_cdfs(mm, ln_pga):
            return get_conditional_fsl_cdfs(
                mm, ln_pga, sigmav, sigmavp, depth, qc, fs, fsl, pa, rd_amax_cap
            )

        fsl_cdfs = pga_grid.quadrature_cdfs(
            conditional_cdfs, m, mu_ln_pga, sigma_ln_pga, truncation_level, n_quadrature
        )
        return fsl_cdfs, ndtri(fsl_cdfs)
    if pga_uncertainty != "closed_form":
        raise ValueError(
            f'pga_uncertainty must be "conditional" or "closed_form", not "{pga_uncertainty}"'
        )
    mu_ln_csr, sigma_ln_csr = get_ln_csr(
        mu_ln_pga, sigma_ln_pga, m, sigmav, sigmavp, depth
    )
    mu_ln_crr, sigma_ln_crr = get_ln_crr(m, mu_ln_csr, qc, fs, sigmavp, pa)
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps = (np.log(fsl) - mu_ln_fsl[:, np.newaxis]) / sigma_ln_fsl[:, np.newaxis]
    fsl_cdfs = ndtr(eps)

    return fsl_cdfs, eps
