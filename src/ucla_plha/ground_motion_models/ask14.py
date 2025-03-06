import numpy as np

def get_im(vs30,rrup,rx,rx1,ry0,m,fault_type,measured_vs30,dip,ztor,**kwargs):
    '''
    Output an array of mean and standard deviation values of the natural log of peak horizontal acceleration.
    vs30 = time-averaged shear wave velocity in upper 30m [m/s]
    rjb = Joyner-Boore source-to-site distance [km]
    rrup = rupture distance [km]
    rx = closest distance to vertical projection of top of fault, measured perpendicular to strike
    rx1 = closest distance to vertical projection of bottom of fault, measured perpendicular to strike
    ry0 = closest distance to rupture, measured parallel to strike
    m = moment magnitude
    fault_type = style of faulting based on rake. 1 = reverse, 2 = normal, 3 = strike slip 
    measured_vs30 = boolean field indicating whether vs30 is measured or inferred
    dip = dip angle in degrees
    ztor = depth to top of rupture [km]

    Keyword arguments (**kwargs)
    z1 = depth to VS = 1.0 k/ms [km]

    Note: We have assumed all earthquakes are mainshocks, and have not included regionalization terms since
    sites are all in California
    '''

    # Define flags
    FRV = np.zeros(len(m), dtype=float)
    FRV[fault_type == 1] = 1.0
    FN = np.zeros(len(m), dtype=float)
    FN[fault_type == 2] = 1.0
    FHW = np.zeros(len(m), dtype=float)
    FHW[(dip != 90) & (rx > rx1)] = 1.0

    # Define z1.0
    z1 = kwargs.get('z1', None)
    c4 = 4.5
    M1 = 6.75
    M2 = 5.0
    a1 = 0.587
    a2 = -0.790
    a3 = 0.275
    a4 = -0.1
    a5 = -0.41
    a6 = 2.154
    a7 = 0.0  # They didn't specify this anywhere in the paper, but it's in the XLS file
    a8 = -0.015
    a10 = 1.735
    a11 = 0.0
    a12 = -0.1
    a13 = 0.6
    a15 = 1.1
    a17 = -0.0072
    vlin = 660
    b = -1.47
    n = 1.5
    c = 2.4
    a43 = 0.1
    a44 = 0.05
    a45 = 0.0
    a46 = -0.05
    
    if(measured_vs30):
        s1 = 0.741
        s2 = 0.501
    else:
        s1 = 0.754
        s2 = 0.520
    s3 = 0.47
    s4 = 0.36

    # Basic form
    f1 = np.empty(len(m), dtype=float)
    f1_filt1 = m >= M1
    f1_filt2 = (M2 <= m) & (m < M1)
    f1_filt3 = m < M2
    c4m = np.empty(len(m))
    m_filt1 = m > 5
    m_filt2 = (4 < m) & (m <= 5)
    m_filt3 = m <= 4
    c4m[m_filt1] = c4
    c4m[m_filt2] = c4 - (c4 - 1) * (5 - m[m_filt2])
    c4m[m_filt3] = 1.0
    R = np.sqrt(rrup ** 2 + c4m ** 2)
    f1[f1_filt1] = a1 + a5 * (m[f1_filt1] - M1) + a8 * (8.5 - m[f1_filt1]) ** 2 + (a2 + a3 * (m[f1_filt1] - M1)) * np.log(R[f1_filt1]) + a17 * rrup[f1_filt1]
    f1[f1_filt2] = a1 + a4 * (m[f1_filt2] - M1) + a8 * (8.5 - m[f1_filt2]) ** 2 + (a2 + a3 * (m[f1_filt2] - M1)) * np.log(R[f1_filt2]) + a17 * rrup[f1_filt2]
    f1[f1_filt3] = a1 + a4 * (M2 - M1) + a8 * (8.5 - M2) ** 2 + a6 * (m[f1_filt3] - M2) + a7 * (m[f1_filt3] - M2) ** 2 + (a2 + a3 * (M2 - M1)) * np.log(R[f1_filt3]) + a17 * rrup[f1_filt3]
    
    # Style of faulting model
    f7 = np.empty(len(m), dtype=float)
    f7[m_filt1] = a11
    f7[m_filt2] = a11 * (m[m_filt2] - 4.0)
    f7[m_filt3] = 0.0

    f8 = np.empty(len(m), dtype=float)
    f8[m_filt1] = a12
    f8[m_filt2] = a12 * (m[m_filt2] - 4.0)
    f8[m_filt3] = 0.0
    
    # Hanging wall model
    T1 = np.zeros(len(m), dtype=float)
    a2HW = 0.2
    r1 = np.abs(rx - rx1)
    r2 = 3 * r1
    h1 = 0.25
    h2 = 1.5
    h3 = -0.75
    ry1 = rx * np.tan(20 * np.pi / 180.0)
    T1[dip>30] = (90 - dip) / 45
    T1[dip<=30] = 60 / 45
    T2 = np.zeros(len(m), dtype=float)
    T2[m >= 6.5] = 1 + a2HW * (m[m >= 6.5] - 6.5)
    T2[(5.5 < m) & (m < 6.5)] = 1 + a2HW * (m[(5.5 < m) & (m < 6.5)] - 6.5) - (1 - a2HW) * (m[(5.5 < m) & (m < 6.5)] - 6.5) ** 2.0
    T2[m >= 5.5] = 0
    T3 = np.zeros(len(m), dtype=float)
    T3[rx < r1] = h1 + h2 * (rx[rx < r1] / r1[rx < r1]) + h3 * (rx[rx < r1] / r1[rx < r1]) **2
    T3[(r1 <= rx) & (rx <= r2) & (r1 != r2)] = 1 - (rx[(r1 <= rx) & (rx <= r2) & (r1 != r2)] - r1[(r1 <= rx) & (rx <= r2) & (r1 != r2)]) / (r2[(r1 <= rx) & (rx <= r2) & (r1 != r2)] - r1[(r1 <= rx) & (rx <= r2) & (r1 != r2)])
    T3[rx > r2] = 0.0
    T4 = np.zeros(len(m), dtype=float)
    T4[ztor <= 10] = 1 - ztor**2 / 100
    T4[ztor > 10] = 0.0
    T5 = np.zeros(len(m), dtype=float)
    T5[0.0 >= ry0 - ry1] = 1.0
    T5[(0.0 < ry0 - ry1) & (ry0 - ry1 < 5)] = 1 - (ry0[(0.0 < ry0 - ry1) & (ry0 - ry1 < 5)] - ry1[(0.0 < ry0 - ry1) & (ry0 - ry1 < 5)]) / 5.0
    f4 = a13 * T1 * T2 * T3 * T4 * T5
    
    # Depth to rupture model
    f6 = np.empty(len(m), dtype=float)
    f6[ztor < 20] = a15 * ztor / 20.0
    f6[ztor >= 20] = a15

    # Soil depth model. Assume z1 = z1ref if it is not specified
    z1ref = 1/1000 * np.exp(-7.67 / 4 * np.log((vs30**4 + 610**4) / (1360**4 + 610**4)))
    if('z1p0' in kwargs):
        z1 = kwargs.get('z1p0')
    else:
        z1 = z1ref
    f10 = np.empty(len(m), dtype=float)
    if(vs30 <= 200):
        f10 = a43 * np.log((z1 + 0.01) / (z1ref + 0.01))
    elif((200 < vs30) & (vs30 <= 300)):
        f10 = np.log((z1 + 0.01) / (z1ref + 0.01))
    elif((300 < vs30) & (vs30 <= 500)):
        f10 = a45 * np.log((z1 + 0.01) / (z1ref + 0.01))
    else:
        f10 = a46 * np.log((z1 + 0.01) / (z1ref + 0.01))
    
    # Aftershock scaling (assume all earthquakes are mainshocks here)
    FAs = 0.0
    f11 = 0.0

    # Regionalization (assume regional term is zero since we are in California)
    regional = 0.0

   
    # Site response model (do this one last because we need to know the ground motion for VS30 = 1180 m/s to compute it)
    f5 = np.zeros(len(m), dtype=float)
    vs30star = np.empty(len(m), dtype=float)
    v1 = 1500
    vs30star[vs30 < v1] = vs30
    vs30star[vs30 >= v1] = v1
    f5[vs30 >= vlin] = (a10 + b * n) * np.log(vs30star[vs30 >= vlin] / vlin)
    
    # Now compute the rock motion amplitude and then do site response
    f5_1180 = (a10 + b * n) * np.log(1180 / vlin)
    sa1180 = np.exp(f1 + f5_1180 + FRV * f7 + FN * f8 + FAs * f11 + FHW * f4 + f6 + f10 + regional)
    f5[vs30 < vlin] = a10 * np.log(vs30star[vs30 < vlin] / vlin) - b * np.log(sa1180[vs30 < vlin] + c) + b * np.log(sa1180[vs30 < vlin] + c * (vs30star[vs30 < vlin] / vlin) ** n)
    
    # Standard deviation
    phi_al = np.empty(len(m), dtype=float)
    phi_al[m < 4.0] = s1
    phi_al[(4 <= m) & (m <= 6)] = s1 + (s2 - s1) / 2.0 * (m[(4 <= m) & (m <= 6)] - 4.0)
    phi_al[m > 6] = s2

    tau_al = np.empty(len(m), dtype=float)
    tau_al[m < 5] = s3
    tau_al[(5 <= m) & (m <= 7)] = s3 + (s4 - s3) / 2.0 * (m[(5 <= m) & (m <= 7)] - 5.0)
    tau_al[(m > 7)] = s4

    phi_amp = 0.4
    phi_b = np.sqrt(phi_al**2 + phi_amp**2)
    dlnamp_dlnsa1180 = np.empty(len(m), dtype=float)
    if(vs30 >= vlin):
        dlnamp_dlnsa1180 = 0.0
    else:
        dlnamp_dlnsa1180 = -b * sa1180 / (sa1180 + c) + b * sa1180 / (sa1180 + c * (vs30 / vlin) ** n)
    phi = np.sqrt(phi_b ** 2 * (1.0 + dlnamp_dlnsa1180) ** 2 + phi_amp**2)
    tau = tau_al * (1 + dlnamp_dlnsa1180)
    mu = f1 + FRV * f7 + FN * f8 + FAs * f11 + f5 + FHW * f4 + f6 + f10 + regional
    sigma = np.sqrt(phi**2 + tau**2)

    return (mu, sigma)
