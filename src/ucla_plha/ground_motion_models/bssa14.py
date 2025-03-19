import numpy as np

def get_im(vs30,rjb,m,fault_type):
    '''
    Output an array of mean and standard deviation values of the natural log of peak horizontal acceleration.
    vs30 = time-averaged shear wave velocity in upper 30m [m/s]
    rjb = Joyner-Boore source-to-site distance [km]
    m = moment magnitude
    fault_type = style of faulting based on rake
    '''
    # vs30
    # PGAr = PGA at reference site condition (eg. rock)
    # T = period in sec
    e0 = 0.4473
    e1 = 0.4856
    e2 = 0.2459
    e3 = 0.4539
    e4 = 1.431
    e5 = 0.05053
    e6 = -0.1662
    mh = 5.5
    c1 = -1.134
    c2 = 0.1917
    c3 = -0.00809
    mref = 4.5
    rref = 1
    h = 4.5
    deltac3 = 0 # no regional adjustment for California
    c = -0.6
    vc = 1500
    vref = 760
    f1 = 0
    f2 = 0.1
    f3 = 0.1
    f4 = -0.15
    f5 = -0.00701
    f6 = -9.900   # Not used since basin term doesn't influence PGA
    f7 = -9.900   # Not used since basin term doesn't influence PGA
    r1 = 110.000
    r2 = 270.000
    deltaphir = 0.1000
    deltaphiv = 0.070
    v1 = 225.0
    v2 = 300.0
    phi1 = 0.695
    phi2 = 0.495
    tau1 = 0.398
    tau2 = 0.348

    rs = np.zeros(len(fault_type), dtype=int)
    ss = np.zeros(len(fault_type), dtype=int)
    ns = np.zeros(len(fault_type), dtype=int)
    rs[fault_type == 1] = 1
    ns[fault_type == 2] = 1
    ss[fault_type == 3] = 1
        
    fe = e1 * ss + e2 * ns + e3 * rs + e4 * (m - mh) + e5 * (m - mh) ** 2
    fe[m > mh] = e1 * ss[m > mh] + e2 * ns[m > mh] + e3 * rs[m > mh] + e6 * (m[m > mh] - mh)
   
    
    # path function
    r = np.sqrt(rjb ** 2 + h ** 2)  # Equation 4
    fpath= (c1 + c2 * (m - mref)) * np.log(r / rref) + (c3 + deltac3) * (r - rref) # Equation 3

    # site term 
    # Compute PGA for rock conditions, which corresponds to VS30 = 760 m/s (i.e., Flin = 0.0)
    pgar= np.exp(fe + fpath)
    
    # Nonlinear site response term Equation 6
    if vs30 < vc:
        lnflin = c * np.log(vs30 / vref)
    else:
        lnflin = c * np.log(vc / vref)

    f2 = f4 * (np.exp(f5 * (min(vs30, vref) - 360)) - np.exp(f5 * (vref - 360)))   # Equation 8
    lnfnl = f1 + f2 * np.log((pgar + f3) / f3)   # Equation 7
    lnfnl[pgar == 0] = 0.0
    lnfnl[(pgar > 0) & (vs30 >= vref)] = 0    
    lnfs= lnflin + lnfnl  # equation 5 (skipping Fdz term since it doesn't influence PGA)
    
    # Standard deviation
    # tau(M) equation 14 in BSSA 14
    tau = np.full(len(fault_type), tau2)
    tau[m<=4.5] = tau1
    tau[(4.5<m) & (m<5.5)] = tau1 + (tau2 - tau1) * (m[(4.5<m) & (m<5.5)] - 4.5)    

    # phi(M) equation 17 in BSSA14
    phi_m = np.empty(len(fault_type), dtype=float)
    phi_m[m <= 4.5] = phi1
    phi_m[(4.5 < m) & (m < 5.5)] = phi1 + (phi2 - phi1) * (m[(4.5 < m) & (m < 5.5)] - 4.5)
    phi_m[m >= 5.5] = phi2

    # phi(M, RJB) equation 16 in BSSA14
    phi_m_rjb = np.empty(len(fault_type), dtype=float)
    phi_m_rjb[rjb <= r1] = phi_m[rjb <= r1]
    phi_m_rjb[(r1 < rjb) & (rjb <= r2)] = phi_m[(r1 < rjb) & (rjb <= r2)] + deltaphir * np.log(rjb[(r1 < rjb) & (rjb <= r2)] / r1) / np.log(r2 / r1)
    phi_m_rjb[rjb > r2] = phi_m[rjb > r2] + deltaphir

    # phi(M, RJB, VS30) equation 15 in BSSA14
    phi = np.empty(len(fault_type), dtype=float)
    phi[vs30 >= v2] = phi_m_rjb[vs30 >= v2]
    phi[(v1 < vs30) & (vs30 < v2)] = phi_m_rjb[(v1 < vs30) & (vs30 < v2)] - deltaphiv * np.log(v2 / vs30) / np.log(v2 / v1)
    phi[vs30 <= v1] = phi_m_rjb[vs30 <= v1] - deltaphiv
    
    # Compute mu and sigma
    mu=(fe + fpath + lnfs)
    sigma= np.sqrt(phi ** 2 + tau ** 2)
    
    return (mu, sigma)
