import numpy as np
import pandas as pd
import scipy as sp
from scipy.stats import norm
from scipy.special import ndtr
import json, jsonschema
from ucla_plha.liquefaction_models import cetin_et_al_2018, idriss_boulanger_2012, ngl_smt_2024
from ucla_plha.ground_motion_models import ask14, bssa14, cb14, cy14
from ucla_plha.geometry import geometry
from importlib_resources import files
import numba

def get_source_data(source_type, source_model, p_xyz, dist_cutoff, gmms):
    if(source_type == 'fault_source_models'):
        # Read files required by all ground motion models
        path = files('ucla_plha').joinpath('source_models/fault_source_models/' + source_model)
        ruptures = pd.read_pickle(str(path.joinpath('ruptures.pkl')), compression='gzip')
        m = ruptures['m'].values
        fault_type = ruptures['style'].values
        ruptures_segments = pd.read_pickle(str(path.joinpath('ruptures_segments.pkl')), compression='gzip')
        segment_index = ruptures_segments['segment_index'].values
        rate = ruptures['rate'].values
        dip = ruptures['dip'].values
        ztor = ruptures['ztor'].values
        zbor = ruptures['zbor'].values

        # Now read files required by ask14, bssa14, cb14, and / or cy14
        # bssa14: rjb
        # ask14: rrup,rx,rx1,ry0,
        # cb14: rjb,rrup,rx
        # cy14: rjb,rrup,rx
        empty_array = np.empty(len(m))
        if(any(gmm in ['bssa14', 'cb14', 'cy14'] for gmm in gmms)):
            tri_segment_id = np.load(str(path.joinpath('tri_segment_id.npy')))
            tri_rjb = np.load(str(path.joinpath('tri_rjb.npy')))
            rjb_all = geometry.point_triangle_distance(tri_rjb, p_xyz, tri_segment_id)
            ruptures_segments['rjb_all'] = rjb_all[segment_index]
            
        if(any(gmm in ['ask14', 'cb14', 'cy14'] for gmm in gmms)):
            rect_segment_id = np.load(str(path.joinpath('rect_segment_id.npy')))
            tri_rrup = np.load(str(path.joinpath('tri_rrup.npy')))
            rect = np.load(str(path.joinpath('rect_rjb.npy')))
            rrup_all = geometry.point_triangle_distance(tri_rrup, p_xyz, tri_segment_id)
            ruptures_segments['rrup_all'] = rrup_all[segment_index]
            Rx_all, Rx1_all, Ry0_all = geometry.get_Rx_Rx1_Ry0(rect, p_xyz, rect_segment_id)
            ruptures_segments['Rx_all'] = Rx_all[segment_index]
            ruptures_segments['Rx1_all'] = Rx1_all[segment_index]
            ruptures_segments['Ry0_all'] = Ry0_all[segment_index]
            ruptures_segments['dip_all'] = dip[segment_index]
            ruptures_segments['ztor_all'] = ztor[segment_index]
            ruptures_segments['zbor_all'] = zbor[segment_index]

        grouped_ruptures_segments = ruptures_segments.groupby('rupture_index')

        if(any(gmm in ['ask14', 'cb14', 'cy14'] for gmm in gmms)):
            rrup = grouped_ruptures_segments['rrup_all'].min().values
            filter = rrup < dist_cutoff
            rx = grouped_ruptures_segments['Rx_all'].min().values
            rx1 = grouped_ruptures_segments['Rx1_all'].min().values
            ry0 = grouped_ruptures_segments['Ry0_all'].min().values
        else:
            rrup = empty_array
            rx = empty_array
            rx1 = empty_array
            ry0 = empty_array

        if(any(gmm in ['bssa14', 'cb14', 'cy14'] for gmm in gmms)):
            rjb = grouped_ruptures_segments['rjb_all'].min().values
            filter = rjb < dist_cutoff
        else:
            rjb = empty_array
        
                
        return(m[filter], fault_type[filter], rate[filter], rjb[filter], rrup[filter], rx[filter], rx1[filter], ry0[filter], dip[filter], ztor[filter], zbor[filter])
    
    elif(source_type == 'point_source_models'):
        path = files('ucla_plha').joinpath('source_models/point_source_models/' + source_model)
        ruptures = pd.read_pickle(str(path.joinpath('ruptures.pkl')), compression='gzip')
        rate = ruptures['rate'].values
        m = ruptures['m'].values
        fault_type = ruptures['style'].values 
        node_index = np.load(str(path.joinpath('node_index.npy')))
        points = np.load(str(path.joinpath('points.npy')))
        dist_temp = np.empty(np.max(node_index)+1, dtype=float)
        for i, ni in enumerate(node_index):
            dist_temp[ni] = np.sqrt((points[i,0] - p_xyz[0])**2 + (points[i,1] - p_xyz[1])**2 + (points[i,2] - p_xyz[2])**2)
        dist = dist_temp[ruptures['node_index'].values]
        filter = dist < dist_cutoff
        dip = np.empty(len(m), dtype=float)

        # using Kaklamanos et al. 2011 guidance fro unknown dip, ztor, and zbor
        # Note fault_type = 1 reverse, 2 normal, 3 strike slip
        dip[fault_type == 1] = 40.0
        dip[fault_type == 2] = 50.0
        dip[fault_type == 3] = 90.0
        w = 10.0 ** (-0.76 + 0.27 * m)
        w[fault_type == 1] = 10.0 ** (-1.61 + 0.41 * m[fault_type == 1])
        w[fault_type == 2] = 10.0 ** (-1.14 + 0.35 * m[fault_type == 2])
        zhyp = 5.63 + 0.68 * m
        zhyp[fault_type == 1] = 11.24 - 0.2 * m[fault_type == 1]
        zhyp[fault_type == 2] = 7.08 + 0.61 * m[fault_type == 2]
        ztor = zhyp - 0.6 * w * np.sin(dip * np.pi / 180.0)
        ztor[ztor < 0.0] = 0.0
        zbor = ztor + w * np.sin(dip * np.pi / 180.0)

        # For point sources, use same distance for Rx, Rx1, Ry0, which will turn off the hanging wall term
        return(m[filter], fault_type[filter], rate[filter], dist[filter], dist[filter], dist[filter], dist[filter], dist[filter], dip[filter], ztor[filter], zbor[filter])

def get_ground_motion_data(gmm, vs30, fault_type, rjb, rrup, rx, rx1, ry0, m, ztor, zbor, dip, z1p0, z2p5, measured_vs30):
    if(gmm == 'bssa14'):
        mu_ln_pga, sigma_ln_pga = bssa14.get_im(vs30, rjb, m, fault_type)
    elif(gmm == 'cb14'):
        mu_ln_pga, sigma_ln_pga = cb14.get_im(vs30,rjb,rrup,rx,m,fault_type,ztor,zbor,dip,z2p5)
    elif(gmm == 'cy14'):
        mu_ln_pga, sigma_ln_pga = cy14.get_im(vs30,rjb,rrup,rx,m,fault_type,measured_vs30,dip,ztor,z1p0)
    elif(gmm == 'ask14'):
        mu_ln_pga, sigma_ln_pga = ask14.get_im(vs30, rrup, rx, rx1, ry0, m, fault_type, dip, ztor, measured_vs30, z1p0)
    else:
        print('incorrect ground motion model')
    return([mu_ln_pga, sigma_ln_pga])

def get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config):
    if(liquefaction_model=='cetin_et_al_2018'):
        c = config['liquefaction_models']['cetin_et_al_2018']
        return cetin_et_al_2018.get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, c['sigmav'], c['sigmavp'], c['vs12'], c['depth'], c['n160'], c['fc'], fsl, c['pa'])
    elif(liquefaction_model=='idriss_boulanger_2012'):
        c = config['liquefaction_models']['idriss_boulanger_2012']
        return idriss_boulanger_2012.get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, c['sigmav'], c['sigmavp'], c['depth'], c['n160'], c['fc'], fsl, c['pa'])
    elif(liquefaction_model=='ngl_smt_2024'):
        c = config['liquefaction_models']['ngl_smt_2024']
        cpt_df = pd.read_csv(c['cpt_data'])
        header = cpt_df.columns.values.tolist()
        if(c['process_cpt']):
            error = None
            # Check that file headings are OK. Otherwise return error.
            # Input options:
            # 1. [depth, qt, fs] plus gamma and gammaw must be specified in config
            # 2. [depth, qt, fs, sigmav, sigmavp] and no gamma and gammaw specified in config
            required_headers = ['depth', 'qt', 'fs']
            if(not all(item in header for item in required_headers)):
                error = "Your CSV file must have these headings: 'depth', 'qt', 'fs'. "
            if(('sigmav' in header) and ('sigmavp' not in header)):
                if(error):
                    error += "<br>If you specify sigmav, you must also specify sigmavp. "
                else:
                    error = "If you specify sigmav, you must also specify sigmavp. "
            if(('sigmavp' in header) and ('sigmav' not in header)):
                if(error):
                    error += "<br>If you specify sigmavp, you must also specify sigmav. "
                else:
                    error = "If you specify sigmavp, you must also specify sigmav. "
            if('sigmav' not in header):
                if(('gamma' not in c.keys()) or ('gammaw' not in c.keys())):
                    if(error):
                        error += '<br>If you do not specify sigmav and sigmavp in cpt_data CSV file, you must specify gamma and gammaw in the config file.'
                    else:
                        error = 'If you do not specify sigmav and sigmavp in cpt_data CSV file, you must specify gamma and gammaw in the config file.'
                sigmav = cpt_df['depth'].values * c['gamma']
                u = (cpt_df['depth'].values - c['dGWT']) * c['gammaw']
                u[u<c['dGWT']] = 0.0
                sigmavp = sigmav - u
            else:
                sigmav = cpt_df['sigmav'].values
                sigmavp = cpt_df['sigmavp'].values
            if(error):
                print("cpt_data ERROR: ", error)
                return
            ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat = ngl_smt_2024.process_cpt(cpt_df['depth'].values, cpt_df['qt'].values, cpt_df['fs'].values, c['dGWT'], sigmav=sigmav, sigmavp=sigmavp)
        else:
            required_headers = ['ztop', 'zbot', 'qc1Ncs', 'Ic', 'sigmav', 'sigmavp', 'Ksat']
            if(not all(item in header for item in required_headers)):
                error = 'If you specify process_cpt = false, the required headers are ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat'
            ztop = cpt_df['ztop'].values
            zbot = cpt_df['zbot'].values
            qc1Ncs = cpt_df['qc1Ncs'].values
            Ic = cpt_df['Ic'].values
            sigmav = cpt_df['sigmav'].values
            sigmavp = cpt_df['sigmavp'].values
            Ksat = cpt_df['Ksat'].values
                                        
        return ngl_smt_2024.get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, m, ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, Ksat, fsl, 101.325)

@numba.jit(nopython=True)           
def get_disagg(liquefaction_hazards, eps, m, rjb, plha_magnitude_bin_edges, plha_distance_bin_edges, plha_epsilon_bin_edges, source_model_weight, ground_motion_model_weight, liquefaction_model_weight):
    N = eps.shape[0]
    M = eps.shape[1]
    plha_magnitude_bin_center = 0.5 * (plha_magnitude_bin_edges[0:-1] + plha_magnitude_bin_edges[1:])
    plha_distance_bin_center = 0.5 * (plha_distance_bin_edges[0:-1] + plha_distance_bin_edges[1:])
    plha_epsilon_bin_center = 0.5 * (plha_epsilon_bin_edges[0:-1] + plha_epsilon_bin_edges[1:])
    plha_disagg = np.zeros((M, len(plha_magnitude_bin_center), len(plha_distance_bin_center), len(plha_epsilon_bin_center)))
    for i in range(N):
        for j in range(len(plha_magnitude_bin_center)):
            for k in range(len(plha_distance_bin_center)):
                for l in range(len(plha_epsilon_bin_center)):
                    plha_disagg[i, j, k, l] += source_model_weight * ground_motion_model_weight * liquefaction_model_weight * np.sum(liquefaction_hazards[i, (m >= plha_magnitude_bin_edges[j]) & (m < plha_magnitude_bin_edges[j+1]) & (rjb >= plha_distance_bin_edges[k]) & (rjb < plha_distance_bin_edges[k+1]) & (eps[i] >= plha_epsilon_bin_edges[l]) & (eps[i] < plha_epsilon_bin_edges[l+1])])
    return plha_disagg

def get_hazard(config_file):
    # Validate config_file against schema. If ngl_smt_2024 liquefaction model is used, the cpt_data file is 
    # validated in the get_liquefaction_hazards function.
    schema = json.load(open(files('ucla_plha').joinpath('ucla_plha_schema.json')))
    config = json.load(open(config_file))
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError as e:
        print("Config File Error:", e.message)
        return
    
    # Read geometry properties
    latitude = config['geometry']['latitude']
    longitude = config['geometry']['longitude']
    elevation = config['geometry']['elevation']
    point = np.asarray([latitude, longitude, elevation])
    dist_cutoff = config['geometry']['dist_cutoff']
    p_xyz = geometry.point_to_xyz(point)

    # Read output properties
    if("psha" in config['output'].keys()):
        pga = np.asarray(config['output']['psha']['pga'], dtype=float)
        output_psha = True
        if("disaggregation" in config['output']['psha'].keys()):
            output_psha_disaggregation = True
            psha_magnitude_bin_edges = np.asarray(config['output']['psha']['disaggregation']['magnitude_bin_edges'], dtype=float)
            psha_distance_bin_edges = np.asarray(config['output']['psha']['disaggregation']['distance_bin_edges'], dtype=float)
            psha_epsilon_bin_edges = np.asarray(config['output']['psha']['disaggregation']['epsilon_bin_edges'], dtype=float)
            psha_magnitude_bin_center = 0.5 * (psha_magnitude_bin_edges[0:-1] + psha_magnitude_bin_edges[1:])
            psha_distance_bin_center = 0.5 * (psha_distance_bin_edges[0:-1] + psha_distance_bin_edges[1:])
            psha_epsilon_bin_center = 0.5 * (psha_epsilon_bin_edges[0:-1] + psha_epsilon_bin_edges[1:])
            psha_disagg = np.zeros((len(pga), len(psha_magnitude_bin_center), len(psha_distance_bin_center), len(psha_epsilon_bin_center)))
        else:
            output_psha_disaggregation = False
        seismic_hazard = np.zeros(len(pga))
    else:
        output_psha = False
        output_psha_disaggregation = False

    if(config['output']['plha']):
        fsl = np.asarray(config['output']['plha']['fsl'], dtype=float)
        output_plha = True
        if("disaggregation" in config['output']['plha'].keys()):
            output_plha_disaggregation = True
            plha_magnitude_bin_edges = np.asarray(config['output']['plha']['disaggregation']['magnitude_bin_edges'], dtype=float)
            plha_distance_bin_edges = np.asarray(config['output']['plha']['disaggregation']['distance_bin_edges'], dtype=float)
            plha_epsilon_bin_edges = np.asarray(config['output']['plha']['disaggregation']['epsilon_bin_edges'], dtype=float)
            plha_magnitude_bin_center = 0.5 * (plha_magnitude_bin_edges[0:-1] + plha_magnitude_bin_edges[1:])
            plha_distance_bin_center = 0.5 * (plha_distance_bin_edges[0:-1] + plha_distance_bin_edges[1:])
            plha_epsilon_bin_center = 0.5 * (plha_epsilon_bin_edges[0:-1] + plha_epsilon_bin_edges[1:])
            plha_disagg = np.zeros((len(pga), len(plha_magnitude_bin_center), len(plha_distance_bin_center), len(plha_epsilon_bin_center)))
        else:
            output_plha_disaggregation = False
        liquefaction_hazard = np.zeros(len(fsl))

    # Loop over all ground motion models to get list of distance types
    gmms = []
    for gmm in config['ground_motion_models'].keys():
        gmms.append(gmm)

    # Loop over source models. We have fault_source_models and point_source_models, so there are two loops
    for source_model in config['source_models'].keys():
       for fault_source_model in config['source_models'][source_model].keys():
            m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = get_source_data(source_model, fault_source_model, p_xyz, dist_cutoff, gmms)
            source_model_weight = config['source_models'][source_model][fault_source_model]['weight']
            # Loop over ground motion models.
            for ground_motion_model in config['ground_motion_models'].keys():
                # retrieve parameters common to all models
                ground_motion_model_weight = config['ground_motion_models'][ground_motion_model]['weight']
                vs30 = config['ground_motion_models'][ground_motion_model]['vs30']
                z1p0 = None
                z2p5 = None
                measured_vs30 = False
                # retrieve model-specific parameters
                if((ground_motion_model == 'ask14') or (ground_motion_model == 'cy14')):
                    if('measured_vs30' in config['ground_motion_models'][ground_motion_model].keys()):
                        measured_vs30 = config['ground_motion_models'][ground_motion_model]['measured_vs30']
                    if('z1p0' in config['ground_motion_models'][ground_motion_model].keys()):
                        z1p0 = config['ground_motion_models'][ground_motion_model]['z1p0']
                if(ground_motion_model == 'cb14'):
                    if('z2p5' in config['ground_motion_models'][ground_motion_model].keys()):
                        z2p5 = config['ground_motion_models'][ground_motion_model]['z2p5']                                                    
                mu_ln_pga, sigma_ln_pga = get_ground_motion_data(ground_motion_model, vs30, fault_type, rjb, rrup, rx, rx1, ry0, m, ztor, zbor, dip, z1p0, z2p5, measured_vs30)
                
                # Compute seismic hazard if requested in config file
                if(output_psha):
                    seismic_hazards = (1 - ndtr((np.log(pga[:, np.newaxis]) - mu_ln_pga) / sigma_ln_pga)) * rate
                    seismic_hazard += source_model_weight * ground_motion_model_weight * np.sum(seismic_hazards, axis=1)
                    # Compute seismic hazard disaggregation if requested in config file
                    if(output_psha_disaggregation):
                        for i in range(len(pga)):
                            eps = (np.log(pga[i]) - mu_ln_pga) / sigma_ln_pga
                            for j in range(len(psha_magnitude_bin_center)):
                                for k in range(len(psha_distance_bin_center)):
                                    for l in range(len(psha_epsilon_bin_center)):
                                        psha_disagg[i, j, k, l] += source_model_weight * ground_motion_model_weight * np.sum(seismic_hazards[i, (m >= psha_magnitude_bin_edges[j]) & (m < psha_magnitude_bin_edges[j+1]) & (rjb >= psha_distance_bin_edges[k]) & (rjb < psha_distance_bin_edges[k+1]) & (eps >= psha_epsilon_bin_edges[l]) & (eps < psha_epsilon_bin_edges[l+1])])
                        

                # Compute liquefaction hazard if requested in config file
                for liquefaction_model in config['liquefaction_models'].keys():
                    liquefaction_model_weight = config['liquefaction_models'][liquefaction_model]['weight']
                    if(output_plha):
                        liquefaction_hazards, eps = get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config)
                        liquefaction_hazards *= rate[:, np.newaxis]
                        liquefaction_hazard += source_model_weight * ground_motion_model_weight * liquefaction_model_weight * np.sum(liquefaction_hazards, axis=0)
                        eps = eps.T
                        liquefaction_hazards = liquefaction_hazards.T
                        # Compute liquefaction hazard disaggregation if requested in config file
                        # if(output_plha_disaggregation):
                        #     for i in range(len(fsl)):
                        #         for j in range(len(plha_magnitude_bin_center)):
                        #             for k in range(len(plha_distance_bin_center)):
                        #                 for l in range(len(plha_epsilon_bin_center)):
                        #                     plha_disagg[i, j, k, l] += source_model_weight * ground_motion_model_weight * liquefaction_model_weight * np.sum(liquefaction_hazards[i, (m >= plha_magnitude_bin_edges[j]) & (m < plha_magnitude_bin_edges[j+1]) & (rjb >= plha_distance_bin_edges[k]) & (rjb < plha_distance_bin_edges[k+1]) & (eps[i] >= plha_epsilon_bin_edges[l]) & (eps[i] < plha_epsilon_bin_edges[l+1])])
                        if(output_plha_disaggregation):
                            plha_disagg = get_disagg(liquefaction_hazards, eps, m, rjb, plha_magnitude_bin_edges, plha_distance_bin_edges, plha_epsilon_bin_edges, source_model_weight, ground_motion_model_weight, liquefaction_model_weight)
    # Now prepare output
    output = {}
    output['input'] = config
    output['output'] = {}   
    if(output_psha):
        if(output_psha_disaggregation):
            for i in range(len(pga)):
                psha_disagg[i] = psha_disagg[i] / seismic_hazard[i] * 100.0
            output['output']['psha'] = {"PGA": pga, "annual_rate_of_exceedance": seismic_hazard, "disaggregation": psha_disagg}
        else:
            output['output']['psha'] = {"PGA": pga, "annual_rate_of_exceedance": seismic_hazard}
    if(output_plha):
        if(output_plha_disaggregation):
            for i in range(len(fsl)):
                plha_disagg[i] = plha_disagg[i] / liquefaction_hazard[i] * 100.0
            output['output']['plha'] = {"FSL": fsl, "annual_rate_of_nonexceedance": liquefaction_hazard, "disaggregation": plha_disagg}
        else:
            output['output']["plha"] = {"FSL": fsl, "annual_rate_of_nonexceedance": liquefaction_hazard}
    return output
