import numpy as np

def get_im(vs30,rjb,rrup,rx,rx1,m,fault_type,ztor,zbor,dip,**kwargs):
    '''
    Output an array of mean and standard deviation values of the natural log of peak horizontal acceleration.
    vs30 = time-averaged shear wave velocity in upper 30m [m/s]
    rjb = Joyner-Boore source-to-site distance [km]
    m = moment magnitude
    fault_type = style of faulting based on rake
    '''
    c = 1.88
    n = 1.18
    h4 = 1.0
    sj = 0.0
    c0 = -4.416
    c1 = 0.984
    c2 = 0.537
    c3 = -1.499
    c4 = -0.496
    c5 = -2.773
    c6 = 0.248
    c7 = 6.768
    c8 = 0.0
    c9 = -0.212
    c10 = 0.720
    c11 = 1.090
    c12 = 2.186
    c13 = 1.420
    c14 = -0.0064
    c15 = -0.202
    c16 = 0.393
    c17 = 0.0977
    c18 = 0.0333
    c19 = 0.00757
    c20 = -0.0055
    k1 = 865.0
    k2 = -1.186
    k3 = 1.839
    a1 = 0.167
    a2 = 0.167
    h1 = 0.241
    h2 = 1.474
    h3 = -0.715
    h5 = -0.337
    h6 = -0.270
    tau1 = 0.409
    tau2 = 0.322
    phi1 = 0.734
    phi2 = 0.492
    phi_ln_af = 0.300
    sigma_m_lt_4p5 = 0.840
    sigma_m_gt_5p5 = 0.588
    if('z2p5' in kwargs):
        z2p5 = kwargs.get('z2p5')
    else:
        z2p5 = np.exp(7.089 - 1.144 * np.log(vs30))  # Default value for California (Equation 33)

    ##############################################################
    # Magnitude term. Equation 2 in CB14
    ##############################################################
    fmag = np.empty(len(m))
    fmag[m <= 4.5] = c0 + c1 * m[m <= 4.5]
    fmag[(4.5 < m) & (m <= 5.5)] = c0 + c1 * m[(4.5 < m) & (m <= 5.5)] + c2 * (m[(4.5 < m) & (m <= 5.5)] - 4.5)
    fmag[(5.5 < m) & (m <= 6.5)] = c0 + c1 * m[(5.5 < m) & (m <= 6.5)] + c2 * (m[(5.5 < m) & (m <= 6.5)] - 4.5) + c3 * (m[(5.5 < m) & (m <= 6.5)] - 5.5)
    fmag[m > 6.5] = c0 + c1 * m[m > 6.5] + c2 * (m[m > 6.5] - 4.5) + c3 * (m[m > 6.5] - 5.5) + c4 * (m[m > 6.5] - 6.5)

    ###############################################################
    # Geometric attenuation term. Equation 3 in CB14
    ###############################################################
    fdis = (c5 + c6 * m) * np.log(np.sqrt(rrup ** 2 + c7 ** 2))

    ###############################################################
    # Style of faulting term. Equations 4, 5, and 6 in CB 14
    ###############################################################
    frv = np.zeros(len(fault_type), dtype=float)
    fnm = np.zeros(len(fault_type), dtype=float)
    frv[fault_type == 1] = 1.0
    fnm[fault_type == 2] = 1.0
    
    fflt_m = np.empty(len(m), dtype = float)
    fflt_m[m <= 4.5] = 0.0
    fflt_m[(4.5 < m) & (m <= 5.5)] = m[(4.5 < m) & (m <= 5.5)] - 4.5
    fflt_m[m > 5.5] = 1.0
    
    fflt_f = c8 * frv + c9 * fnm
    fflt = fflt_f * fflt_m

    ###############################################################
    # Hanging wall term
    ###############################################################
    r1 = (zbor - ztor) * np.cos(dip * np.pi / 180.0)
    r2 = 62.0 * m - 350.0    # Equation 12
    
    # Equation 16
    f_hng_delta = (90 - dip) / 45

    # Equation 9
    f1_rx = h1 + h2 * (rx / r1) + h3 * (rx / r1) ** 2

    # Equation 10
    f2_rx = h4 + h5 * (rx - r1) / (r2 - r1) + h6 * ((rx - r1)/(r2 - r1)) ** 2

    # Equation 8
    fhng_rx = np.empty(len(m), dtype=float)
    fhng_rx[rx < 0] = 0.0
    fhng_rx[(0 <= rx) & (rx < r1)] = f1_rx[(0 <= rx) & (rx < r1)]
    fhng_rx[rx >= r1] = f2_rx[rx >= r1]
    fhng_rx[(rx >= r1) & (fhng_rx < 0.0)] = 0.0

    # Equation 13
    fhng_rrup = np.empty(len(m), dtype=float)
    fhng_rrup[rrup == 0] = 1.0
    fhng_rrup[rrup > 0] = (rrup[rrup > 0] - rjb[rrup > 0]) / rrup[rrup > 0]

    # Equation 14
    fhng_m = np.empty(len(m), dtype=float)
    fhng_m[m <= 5.5] = 0.0
    fhng_m[(5.5 <= m) & (m <= 6.5)] = (m[(5.5 <= m) & (m <= 6.5)] - 5.5) * (1.0 + a2 * (m[(5.5 <= m) & (m <= 6.5)] - 6.5))
    fhng_m[m > 6.5] = 1.0 + a2 * (m[m > 6.5] - 6.5)    

    # Equation 15
    fhng_z = np.empty(len(m), dtype=float)
    fhng_z[ztor <= 16.66] = 1.0 - 0.06 * ztor[ztor <= 16.66]
    fhng_z[ztor > 16.66] = 0.0

    # Equation 16
    fhng_delta = (90 - dip) / 45.0
    
    # Equation 7
    fhng = c10 * fhng_rx * fhng_rrup * fhng_m * fhng_z * fhng_delta
    
    ####################################################################
    # Basin Response Term
    ####################################################################
    
    # Equation 20
    fsed = np.empty(len(m), dtype=float)
    fsed[z2p5 <= 1.0] = (c14 + c15 * sj) * (z2p5 - 1.0)
    fsed[(1.0 < z2p5) & (z2p5 <= 3.0)] = 0.0
    fsed[z2p5 > 3.0] = c16 * k3 * np.exp(-0.75) * (1.0 - np.exp(-0.25 * (z2p5 - 3.0)))

    ####################################################################
    # Hypocentral Depth Term
    ####################################################################

    # Equation 36
    fdz_m = np.empty(len(m), dtype=float)
    fdz_m[m < 6.75] = -4.317 + 0.984 * m[m < 6.75]
    fdz_m[m >= 6.75] = 2.325

    # Equation 37
    fdz_delta = np.zeros(len(m), dtype=float)
    fdz_delta[dip <= 40.0] = 0.0445 * (dip[dip <= 40.0] - 40.0)

    # Equation 35
    lndz = fdz_m + fdz_delta
    lndz[lndz > np.log(0.9 * (zbor - ztor))] = np.log(0.9 * (zbor[lndz > np.log(0.9 * (zbor - ztor))] - ztor[lndz > np.log(0.9 * (zbor - ztor))]))
    zhyp = ztor + np.exp(lndz)
    
    # Equation 22
    fhyp_h = np.empty(len(m), dtype=float)
    fhyp_h[zhyp <= 7.0] = 0.0
    fhyp_h[(7.0 < zhyp) & (zhyp <= 20.0)] = zhyp - 7.0
    fhyp_h[zhyp > 20] = 13.0

    # Equation 23
    fhyp_m = np.empty(len(m), dtype=float)
    fhyp_m[m <= 5.5] = c17
    fhyp_m[(5.5 < m) & (m <= 6.5)] = c17 + (c18 - c17) * (m[(5.5 < m) & (m <= 6.5)] - 5.5)
    fhyp_m[m > 6.5] = c18

    # Equation 21
    fhyp = fhyp_h * fhyp_m

    ####################################################################
    # Fault Dip Term
    ####################################################################

    # Equation 24
    fdip = np.empty(len(m), dtype=float)
    fdip[m <= 4.5] = c19 * dip[m <= 4.5]
    fdip[(4.5 <= m) & (m <= 5.5)] = c19 * (5.5 - m[(4.5 <= m) & (m <= 5.5)]) * dip[(4.5 <= m) & (m <= 5.5)]
    fdip[m > 5.5] = 0.0

    ####################################################################
    # Anelastic Attenuation Term
    ####################################################################

    # Equation 25
    fatn = np.zeros(len(m), dtype=float)
    fatn[rrup > 80] = c20 * (rrup[rrup > 80] - 80)

    ####################################################################
    # Shallow Site Response Term (do this last because we need A1100)
    ####################################################################
    
    # A1100 is the value of PGA for VS30 = 1100 m/s, which can be obtained by summing all other terms with site term = 1.0
    a1100 = np.exp(fmag + fdis + fflt + fhng + fsed + fhyp + fdip + fatn)

    # Equation 18
    if(vs30 <= k1):
        f_site_g = c11 * np.log(vs30 / k1) + k2 * (np.log(a1100 + c * (vs30 / k1) ** n) - np.log(a1100 + c))
    else:
        f_site_g = (c11 + k2 * n) * np.log(vs30 / k1)

    # Equation 17 (No Japan term here because we are in California)
    fsite = f_site_g  

    # Equation 1
    mu = fmag + fdis + fflt + fhng + fsite + fsed + fhyp + fdip + fatn

    #####################################################################
    # Aleatory Variability Model
    #####################################################################

    # Equation 27
    tau_lny = np.empty(len(m), dtype=float)
    tau_lny[m <= 4.5] = tau1
    tau_lny[(4.5 < m) & (m < 5.5)] = tau2 + (tau1 - tau2) * (5.5 - m[(4.5 < m) & (m < 5.5)])
    tau_lny[m >= 5.5] = tau2

    # Equation 28
    phi_lny = np.empty(len(m), dtype=float)
    phi_lny[m <= 4.5] = phi1
    phi_lny[(4.5 < m) & (m < 5.5)] = phi2 + (phi1 - phi2) * (5.5 - m[(4.5 < m) & (m < 5.5)])
    phi_lny[m >= 5.5] = phi2

    # Equation 31
    if(vs30 >= k1):
        alpha = np.zeros(len(m), dtype=float)
    else:
        alpha = k2 * a1100 * (a1100 + c * (vs30 / k1) ** n) ** -1 - (a1100 + c) ** -1

    # Equation 29
    tau_lnyb = tau_lny   # Our IM is PGA, so these are the same
    tau_lnpgab = tau_lny
    phi_lnyb = np.sqrt(phi_lny ** 2 + phi_ln_af **2)
    phi_lnpgab = np.sqrt(phi_lny**2 - phi_ln_af**2)
    tau = np.sqrt(tau_lnyb ** 2 + alpha ** 2 * tau_lnpgab ** 2 + 2 * alpha * tau_lnyb * tau_lnpgab)
    phi = np.sqrt(phi_lnyb ** 2 + phi_ln_af**2 + alpha ** 2 * phi_lnyb ** 2 + 2 * alpha * phi_lnyb * phi_lnpgab)

    # Equation 32
    sigma = np.sqrt(tau ** 2 + phi ** 2)

    return (mu, sigma)
