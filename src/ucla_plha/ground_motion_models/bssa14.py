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
    e0=0.4473
    e1=0.4856
    e2=0.2459
    e3=0.4539
    e4=1.431
    e5=0.05053
    e6=-0.1662
    mh= 5.5
    c1=-1.134
    c2=0.1917
    c3=-0.00809
    mref=4.5
    rref=1
    h=4.5
    deltac3=0 #check for california, maybe nonzero
    c=-0.6
    vc=1500
    vref=760
    f1=0
    f2=0.1
    f3 = 0.1
    f4=-0.15
    f5=-0.00701
    f6=-9.900
    f7=-9.900
    r1=110.000
    r2=270.000
    deltaphir=0.1000
    deltaphiv=.070
    v1=225
    v2=300
    phi1=0.695
    phi2=0.495
    tau1=0.398
    tau2=0.348

    rs = np.zeros(len(fault_type), dtype=int)
    ss = np.zeros(len(fault_type), dtype=int)
    ns = np.zeros(len(fault_type), dtype=int)
    rs[fault_type==1] = 1
    ns[fault_type==2] = 1
    ss[fault_type==3] = 1
        
    fe = e1*ss+e2*ns+e3*rs+e4*(m-mh)+e5*(m-mh)**2
    fe[m>mh] = e1*ss[m>mh]+e2*ns[m>mh]+e3*rs[m>mh]+e6*(m[m>mh]-mh)
   
    
    # path function
    r=np.sqrt(rjb**2+h**2)

    fpath= (c1+c2*(m-mref))*np.log(r/rref)+(c3+deltac3)*(r-rref)

    # site term 
    pgar= np.exp(fe+fpath)
    
    if vs30 < vc:
        lnflin = c * np.log(vs30 / vref)
    else:
        lnflin = c * np.log(vc / vref)

    # Nonlinear Site Term
    f2 = f4 * (np.exp(f5 * (min(vs30, vref) - 360)) - np.exp(f5 * (vref - 360)))
    lnfnl = f1+f2 * np.log((pgar + f3) / f3)
    lnfnl[pgar == 0] = 0.0
    lnfnl[(pgar > 0) & (vs30 >= vref)] = 0    
    lnfs= lnflin + lnfnl
    
    # standard deviation, sigma

    tau = np.full(len(fault_type), tau2)
    tau[m<=4.5] = tau1
    tau[(4.5<m) & (m<5.5)] = tau1 + (tau2 - tau1) * (m[(4.5<m) & (m<5.5)] - 4.5)
        
    phi = np.empty(len(fault_type))
    cond1 = (vs30 >= v2) & (rjb <= r1) & (m <= 4.5)
    phi[cond1] = phi1
    cond2 = (vs30 >= v2) & (rjb <= r1) & (4.5 < m) & (m < 5.5)
    phi[cond2] = phi1 + (phi2 - phi1) * (m[cond2] - 4.5)
    cond3 = (vs30 >= v2) & (rjb <= r1) & (m >= 5.5)
    phi[cond3] = phi2
    cond4 = (vs30 >= v2) & (r1 < rjb) & (rjb <= r2) & (m <= 4.5)
    phi[cond4] = phi1 + deltaphir * (np.log(rjb[cond4]/r1)/np.log(r2/r1))
    cond5 = (vs30 >= v2) & (r1 < rjb) & (rjb <= r2) & (4.5 < m) & (m < 5.5)
    phi[cond5] = phi1+(phi2-phi1)*(m[cond5]-4.5)+deltaphir*(np.log(rjb[cond5]/r1)/np.log(r2/r1))
    cond6 = (vs30 >= v2) & (r1 < rjb) & (rjb <= r2) & (m >= 5.5)
    phi[cond6] = phi2+deltaphir*(np.log(rjb[cond6]/r1)/np.log(r2/r1))
    cond7 = (vs30 >= v2) & (rjb > r2) & (m <= 4.5)
    phi[cond7] = phi1+deltaphir
    cond8 = (vs30 >= v2) & (rjb > r2) & (4.5 < m) & (m < 5.5)
    phi[cond8] = phi1+(phi2-phi1)*(m[cond8]-4.5)+deltaphir
    cond9 = (vs30 >= v2) & (rjb > r2) & (m >= 5.5)
    phi[cond9] = phi2+deltaphir

    cond10 = (v1 <= vs30) & (vs30 <= v2) & (rjb <= r1) & (m <= 4.5)
    phi[cond10]= phi1-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond11 = (v1 <= vs30) & (vs30 <= v2) & (rjb <= r1) & (4.5 < m) & (m < 5.5)
    phi[cond11] = phi1+(phi2-phi1)*(m[cond11]-4.5)-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond12 = (v1 <= vs30) & (vs30 <= v2) & (rjb <= r1) & (m > 5.5)
    phi[cond12] = phi2-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond13 = (v1 <= vs30) & (vs30 <= v2) & (r1 < rjb) & (rjb <= r2) & (m <= 4.5)
    phi[cond13] = phi1+deltaphir*(np.log(rjb[cond13]/r1)/np.log(r2/r1))-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond14 = (v1 <= vs30) & (vs30 <= v2) & (r1 < rjb) & (rjb <= r2) & (4.5 < m) & (m < 5.5)
    phi[cond14] = phi1+(phi2-phi1)*(m[cond14]-4.5)+deltaphir*(np.log(rjb[cond14]/r1)/np.log(r2/r1))-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond15 = (v1 <= vs30) & (vs30 <= v2) & (r1 < rjb) & (rjb <= r2) & (m > 5.5)
    phi[cond15] = phi2+deltaphir*(np.log(rjb[cond15]/r1)/np.log(r2/r1))-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond16 = (v1 <= vs30) & (vs30 <= v2) & (rjb > r2) & (m <= 4.5)
    phi[cond16] = phi1+deltaphir-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond17 = (v1 <= vs30) & (vs30 <= v2) & (rjb > r2) & (m > 5.5)
    phi[cond17] = phi1+(phi2-phi1)*(m[cond17]-4.5)+deltaphir-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))
    cond18 = (v1 <= vs30) & (vs30 <= v2) & (rjb > r2) & (m > 5.5)
    phi[cond18] = phi2+deltaphir-deltaphiv*(np.log(v2/vs30)/np.log(v2/v1))

    cond19 = (vs30 < v1) & (rjb <= r1) & (m <= 4.5)
    phi[cond19] = phi1-deltaphiv
    cond20 = (vs30 < v1) & (rjb <= r1) & (4.5 < m) & (m < 5.5)
    phi[cond20] = phi1+(phi2-phi1)*(m[cond20]-4.5)-deltaphiv
    cond21 = (vs30 < v1) & (rjb <= r1) & (m >= 5.5)
    phi[cond21] = phi2-deltaphiv
    cond22 = (vs30 < v1) & (r1 < rjb) & (rjb <= r2) & (m <= 4.5)
    phi[cond22] = phi1+deltaphir*(np.log(rjb[cond22]/r1)/np.log(r2/r1))-deltaphiv
    cond23 = (vs30 < v1) & (r1 < rjb) & (rjb <= r2) & (4.5 < m) & (m < 5.5)
    phi[cond23] = phi1+(phi2-phi1)*(m[cond23]-4.5)+deltaphir*(np.log(rjb[cond23]/r1)/np.log(r2/r1))-deltaphiv
    cond24 = (vs30 < v1) & (r1 < rjb) & (rjb <= r2) & (m >= 5.5)
    phi[cond24] = phi2+deltaphir*(np.log(rjb[cond24]/r1)/np.log(r2/r1))-deltaphiv
    cond25 = (vs30 < v1) & (rjb > r2) & (m <= 4.5)
    phi[cond25] = phi1+deltaphir-deltaphiv
    cond26 = (vs30 < v1) & (rjb > r2) & (4.5 < m) & (m < 5.5)
    phi[cond26] = phi1+(phi2-phi1)*(m[cond26]-4.5)+deltaphir-deltaphiv
    cond27 = (vs30 < v1) & (rjb > r2) & (m >= 5.5)
    phi[cond27] = phi2+deltaphir-deltaphiv
    
    sigma= np.sqrt(phi**2+tau**2)

    mu=(fe+fpath+lnfs)
    
    return (mu, sigma)
