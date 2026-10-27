import numpy as np
from scipy.stats import norm
from scipy.special import ndtr

from ucla_plha import pga_grid


def pdf(x, mu, sigma):
    return 1.0 / np.sqrt(2.0 * np.pi) / sigma * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def get_ln_crr(mu_ln_pga, m, qc1Ncs, sigmavp, pa):
    """
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
    """
    lambda_dr = 1.798
    dr = 47.8 * qc1Ncs**0.264 - 106.3
    dr[dr < 0] = 0
    dr[dr > 160] = 160
    dr_hat = (dr**lambda_dr - 1) / lambda_dr
    ksigma = (sigmavp / pa) ** (-0.0015 * dr)
    ksigma[ksigma > 1.2] = 1.2
    neq = np.exp(0.46305 - 0.4082 * np.exp(mu_ln_pga) + 0.2332 * m)
    msf = (14 / neq) ** 0.2
    mu_ln_crr = (
        -2.799
        + 0.00092 * dr_hat[np.newaxis, :]
        + np.log(ksigma[np.newaxis, :])
        + np.log(msf[:, np.newaxis])
    )
    sigma_ln_crr = 0.44

    return mu_ln_crr, sigma_ln_crr


def get_ln_csr(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, sigmav, sigmavp):
    """
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
    """
    d = (ztop + zbot) / 2.0
    alpha = np.exp(-4.373 + 0.4491 * m)
    beta = -20.11 + 6.247 * m
    rd = (1.0 - alpha[:, np.newaxis]) * np.exp(-d / beta[:, np.newaxis]) + alpha[
        :, np.newaxis
    ]
    mu_ln_csr = mu_ln_pga[:, np.newaxis] + np.log(0.65 * sigmav / sigmavp * rd)
    sigma_ln_csr = sigma_ln_pga

    return mu_ln_csr, sigma_ln_csr


def get_fsl_cdfs(
    mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa,
    return_eps=True, pga_uncertainty="shared", truncation_level=None, n_quadrature=128,
):
    """
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
    return_eps = if False, eps is not computed (it is used only for disaggregation) and None is
        returned for it
    pga_uncertainty = treatment of the uncertainty of the peak acceleration, string:
        "shared" (default): all layers of the profile see the same peak acceleration. The
            probability of profile manifestation is computed conditioned on PGA (with the
            magnitude scaling factor of the conditioned PGA) and integrated over the distribution
            of ln(PGA) of the event (Gauss-Legendre quadrature with n_quadrature points):
            P(FS < fsl | event) = integral of [1 - prod_j (1 - pfmt_j pfts_j(m, PGA) pfs_j Ksat_j)^(t_j / tc)]
            f(ln PGA | mu_ln_pga, sigma_ln_pga) d ln PGA, with pfts_j from sigma_ln_crr = 0.44 only.
        "per_layer": ucla_plha 2.x and earlier. The PGA uncertainty is included in each layer
            separately (sigma_ln_fsl = sqrt(0.44^2 + sigma_ln_pga^2), magnitude scaling factor of
            the median PGA) before the product over the layers, as if each layer had an
            independent PGA. truncation_level is not used.
    truncation_level = upper truncation of the distribution of ln(PGA) in standard deviations
        ("shared" only), as in the PSHA (plha.get_exceedance); None (default) is untruncated
    n_quadrature = number of quadrature points ("shared")

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
    """
    if pga_uncertainty == "shared":
        fsl_cdfs = _get_fsl_cdfs_quadrature(
            mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa,
            truncation_level, n_quadrature,
        )
        return fsl_cdfs, (get_eps(fsl_cdfs) if return_eps else None)
    if pga_uncertainty != "per_layer":
        raise ValueError(
            f'pga_uncertainty must be "shared" or "per_layer", not "{pga_uncertainty}"'
        )
    tc = 2.0
    mu_ln_crr, sigma_ln_crr = get_ln_crr(mu_ln_pga, m, qc1Ncs, sigmavp, pa)
    mu_ln_csr, sigma_ln_csr = get_ln_csr(
        mu_ln_pga, sigma_ln_pga, m, ztop, zbot, sigmav, sigmavp
    )
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    sigma_ln_fsl = np.sqrt(sigma_ln_crr**2 + sigma_ln_csr**2)
    eps_triggering = (
        np.log(fsl[np.newaxis, np.newaxis, :]) - mu_ln_fsl[:, :, np.newaxis]
    ) / sigma_ln_fsl[:, np.newaxis, np.newaxis]
    pfts = ndtr(eps_triggering)
    pfs = 1.0 - 1.0 / (1.0 + np.exp(-14.8 * (Ic / 2.635 - 1.0)))
    pfmt = 1.0 / (1.0 + np.exp(-(7.613 - 0.338 * ztop - 3.042 * Ic)))
    t = zbot - ztop
    fsl_cdfs = 1.0 - np.prod(
        (
            1.0
            - pfmt[np.newaxis, :, np.newaxis]
            * pfts
            * pfs[np.newaxis, :, np.newaxis]
            * Ksat[np.newaxis, :, np.newaxis]
        )
        ** (t[np.newaxis, :, np.newaxis] / tc),
        axis=1,
    )
    if not return_eps:
        return fsl_cdfs, None
    eps_manifestation = norm.ppf(fsl_cdfs)

    return fsl_cdfs, eps_manifestation


def get_eps(fsl_cdfs):
    """Epsilon of profile manifestation (see get_fsl_cdfs): the inverse standard normal CDF of
    the cumulative probabilities, same shape."""
    return norm.ppf(fsl_cdfs)


def get_conditional_fsl_cdfs(m, ln_pga, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa):
    """
    Cumulative distribution functions of the factor of safety against profile manifestation
    conditioned on magnitude and peak acceleration (the same PGA in all layers).

    Inputs:
    m = magnitude, Numpy array, dtype=float, length = K
    ln_pga = natural log of peak acceleration [g], Numpy array, dtype=float, length = K
    ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa: as for get_fsl_cdfs

    Outputs:
    fsl_cdfs = P(FS < fsl | m, PGA) = 1 - prod_j (1 - pfmt_j pfts_j pfs_j Ksat_j)^(t_j / tc), with
        pfts_j = Phi((ln fsl - mu_ln_crr_j(m, PGA) + ln CSR_j(m, PGA)) / 0.44), Numpy ndarray,
        dtype=float, shape = K x L

    Notes:
    The magnitude scaling factor is that of the given PGA. Layers with pfmt pfs Ksat = 0 (e.g.
    above the groundwater table) have no effect and are skipped.
    """
    tc = 2.0
    m = np.asarray(m, dtype=float)
    ln_pga = np.asarray(ln_pga, dtype=float)
    pfs = 1.0 - 1.0 / (1.0 + np.exp(-14.8 * (Ic / 2.635 - 1.0)))
    pfmt = 1.0 / (1.0 + np.exp(-(7.613 - 0.338 * ztop - 3.042 * Ic)))
    a = pfmt * pfs * Ksat
    keep = a != 0.0
    fsl = np.asarray(fsl, dtype=float)
    if not np.any(keep):
        return np.zeros((len(m), len(fsl)))
    t = (zbot - ztop)[keep]
    mu_ln_crr, sigma_ln_crr = get_ln_crr(ln_pga, m, qc1Ncs[keep], sigmavp[keep], pa)
    mu_ln_csr, _ = get_ln_csr(
        ln_pga, 0.0, m, ztop[keep], zbot[keep], sigmav[keep], sigmavp[keep]
    )
    mu_ln_fsl = mu_ln_crr - mu_ln_csr
    pfts = ndtr(
        (np.log(fsl)[np.newaxis, np.newaxis, :] - mu_ln_fsl[:, :, np.newaxis]) / sigma_ln_crr
    )
    log_survival = np.sum(
        (t / tc)[np.newaxis, :, np.newaxis]
        * np.log1p(-a[keep][np.newaxis, :, np.newaxis] * pfts),
        axis=1,
    )
    return -np.expm1(log_survival)


def _get_fsl_cdfs_quadrature(
    mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa,
    truncation_level=None, n_quadrature=128,
):
    """P(FS < fsl | event) with the PGA shared by the layers, by Gauss-Legendre quadrature of
    the conditional probabilities over the (truncated) normal distribution of ln(PGA) of each
    event (pga_grid.quadrature_cdfs). Shape N x L."""

    def conditional_cdfs(mm, ln_pga):
        return get_conditional_fsl_cdfs(
            mm, ln_pga, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, pa
        )

    return pga_grid.quadrature_cdfs(
        conditional_cdfs, m, mu_ln_pga, sigma_ln_pga, truncation_level, n_quadrature
    )
