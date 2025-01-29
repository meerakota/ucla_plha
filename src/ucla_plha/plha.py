import numpy as np
import pandas as pd
import scipy as sp
from scipy.special import ndtr
import json, jsonschema
from ucla_plha.liquefaction_models import cetin_et_al_2018, idriss_boulanger_2012, ngl_smt_2024
from ucla_plha.ground_motion_models import ask14, bssa14, cb14, cy14
from ucla_plha.geometry import geometry
from importlib_resources import files

def get_source_data(source_type, source_model, dist_type, p_xyz, dist_cutoff):
    if(source_type == 'fault_source_models'):
        path = files('ucla_plha').joinpath('source_models/fault_source_models/' + source_model)
        ruptures = pd.read_pickle(str(path.joinpath('ruptures.pkl')), compression='gzip')
        fault_id = np.load(str(path.joinpath('fault_id.npy')))
        m = ruptures['m'].values
        fault_type = ruptures['style'].values
        ruptures_segments = pd.read_pickle(str(path.joinpath('ruptures_segments.pkl')), compression='gzip')
        segment_index = ruptures_segments['segment_index'].values
        if(dist_type == 'rjb'):
            tri = np.load(str(path.joinpath('tri_rjb.npy')))
        elif(dist_type == 'rrup'):
            tri = np.load(str(path.joinpath('tri_rrup.npy')))
        dist_all = geometry.point_triangle_distance(tri, p_xyz, fault_id)
        ruptures_segments['dist_all'] = dist_all[segment_index]
        dist = ruptures_segments.groupby('rupture_index')['dist_all'].min().values
        filter = (dist < dist_cutoff)
        rate = ruptures['rate'].values
    elif(source_type == 'point_source_models'):
        path = files('ucla_plha').joinpath('source_models/point_source_models/' + source_model)
        ruptures = pd.read_pickle(str(path.joinpath('ruptures.pkl')), compression='gzip')
        m = ruptures['m'].values
        fault_type = ruptures['style'].values 
        node_index = np.load(str(path.joinpath('node_index.npy')))
        points = np.load(str(path.joinpath('points.npy')))
        dist_temp = np.empty(np.max(node_index)+1, dtype=float)
        for i, ni in enumerate(node_index):
            dist_temp[ni] = np.sqrt((points[i,0] - p_xyz[0])**2 + (points[i,1] - p_xyz[1])**2 + (points[i,2] - p_xyz[2])**2)
        dist = dist_temp[ruptures['node_index'].values]
        filter = (dist < dist_cutoff)
        rate = ruptures['rate'].values
    return([m[filter], fault_type[filter], rate[filter], dist[filter]])

def get_ground_motion_data(gmm, vs30, dist, m, fault_type):
    if(gmm == 'bssa14'):
        mu_ln_pga, sigma_ln_pga = bssa14.get_im(vs30, dist, m, fault_type)
    elif(gmm == 'cb14'):
        mu_ln_pga, sigma_ln_pga = cb14.get_im(vs30, dist, m, fault_type)
    elif(gmm == 'ba14'):
        mu_ln_pga, sigma_ln_pga = cy14.get_im(vs30, dist, m, fault_type)
    elif(gmm == 'ask14'):
        mu_ln_pga, sigma_ln_pga = ask14.get_im(vs30, dist, m, fault_type)
    else:
        print('incorrect ground motion model')
    return([mu_ln_pga, sigma_ln_pga])

def get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config):
    if(liquefaction_model=='cetin_et_al_2018'):
        c = config['liquefaction_models']['cetin_et_al_2018']
        return cetin_et_al_2018.get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, fsl, m, c['sigmav'], c['sigmavp'], c['vs12'], c['depth'], c['n160'], c['fc'], c['pa'])
    elif(liquefaction_model=='idriss_boulanger_2012'):
        c = config['liquefaction_models']['idriss_boulanger_2012']
        return idriss_boulanger_2012.get_fsl_cdfs(mu_ln_pga, sigma_ln_pga, fsl, m, c['sigmav'], c['sigmavp'], c['depth'], c['n160'], c['fc'], c['pa'])
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
            ztop, zbot, qc1Ncs_lay, Ic_lay, sigmav_lay, sigmavp_lay, Ksat_lay = ngl_smt_2024.process_cpt(cpt_df['depth'].values, cpt_df['qt'].values, cpt_df['fs'].values, c['dGWT'], sigmav=sigmav, sigmavp=sigmavp)
        else:
            required_headers = ['ztop', 'zbot', 'qc1Ncs_lay', 'Ic_lay', 'sigmav_lay', 'sigmavp_lay', 'Ksat_lay']
            if(not all(item in header for item in required_headers)):
                error = 'If you specify process_cpt = false, the required headers are ztop, zbot, qc1Ncs_lay, Ic_lay, sigmav_lay, sigmavp_lay, Ksat_lay'
            ztop = cpt_df['ztop'].values
            zbot = cpt_df['zbot'].values
            qc1Ncs_lay = cpt_df['qc1Ncs_lay'].values
            Ic_lay = cpt_df['Ic_lay'].values
            sigmav_lay = cpt_df['sigmav_lay'].values
            sigmavp_lay = cpt_df['sigmavp_lay'].values
            Ksat_lay = cpt_df['Ksat_lay'].values

        return ngl_smt_2024.get_fsl_cdfs(ztop, zbot, qc1Ncs_lay, Ic_lay, sigmav_lay, sigmavp_lay, Ksat_lay, mu_ln_pga, sigma_ln_pga, m, fsl, N = c['N'])
            

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

    dist_types = {'bssa14': 'rjb', 'cy14': 'rjb', 'cb14': 'rjb', 'as14': 'rjb'}
    
    # Read geometry properties
    latitude = config['geometry']['latitude']
    longitude = config['geometry']['longitude']
    elevation = config['geometry']['elevation']
    point = np.asarray([latitude, longitude, elevation])
    rjb_cutoff = config['geometry']['rjb_cutoff']
    rrup_cutoff = config['geometry']['rrup_cutoff']
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
            output_psha_disaggregation = True
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
    distance_types = []
    for gmm in config['ground_motion_models'].keys():
        if(dist_types[gmm] not in distance_types):
            distance_types.append(dist_types[gmm])

    # Loop over source models. This is nested between fault_source_models and point_source_models, so there are two loops
    for source_model in config['source_models'].keys():
        for fault_source_model in config['source_models'][source_model].keys():
            source_model_weight = config['source_models'][source_model][fault_source_model]['weight']
            # Compute distances we need for the ground motion models we are using
            if('rjb' in distance_types):
                m, fault_type, rate, rjb = get_source_data(source_model, fault_source_model, 'rjb', p_xyz, rjb_cutoff)
            if('rrup' in distance_types):
                m, fault_type, rate, rrup = get_source_data(source_model, fault_source_model, 'rrup', p_xyz, rrup_cutoff)
            # Loop over ground motion models.
            for ground_motion_model in config['ground_motion_models'].keys():
                ground_motion_model_weight = config['ground_motion_models'][ground_motion_model]['weight']
                vs30 = config['ground_motion_models'][ground_motion_model]['vs30']
                if(dist_types[ground_motion_model] == 'rjb'):
                    mu_ln_pga, sigma_ln_pga = get_ground_motion_data(ground_motion_model, vs30, rjb, m, fault_type)
                elif(dist_types[ground_motion_model] == 'rrup'):
                    mu_ln_pga, sigma_ln_pga = get_ground_motion_data(ground_motion_model, vs30, rrup, m, fault_type)
                else:
                    print('incorrect ground motion model')
                # Compute seismic hazard if requested in config file
                if(output_psha):
                    seismic_hazards = (1 - ndtr((np.log(pga[:, np.newaxis]) - mu_ln_pga) / sigma_ln_pga)) * rate
                    seismic_hazard += source_model_weight * ground_motion_model_weight * np.sum(seismic_hazards, axis=1)
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
                        if(output_plha_disaggregation):
                            for i in range(len(fsl)):
                                for j in range(len(plha_magnitude_bin_center)):
                                    for k in range(len(plha_distance_bin_center)):
                                        for l in range(len(plha_epsilon_bin_center)):
                                            plha_disagg[i, j, k, l] += source_model_weight * ground_motion_model_weight * np.sum(seismic_hazards[i, (m >= plha_magnitude_bin_edges[j]) & (m < plha_magnitude_bin_edges[j+1]) & (rjb >= plha_distance_bin_edges[k]) & (rjb < plha_distance_bin_edges[k+1]) & (eps[i] >= plha_epsilon_bin_edges[l]) & (eps[i] < plha_epsilon_bin_edges[l+1])])

    # Now prepare output
    output = {}
    output['input'] = config
    output['output'] = {}   
    if(output_psha):
        if(output_psha_disaggregation):
            for i in range(len(pga)):
                psha_disagg[i] = psha_disagg[i] / seismic_hazard[i] * 100
            output['output']['psha'] = {"PGA": pga, "annual_rate_of_exceedance": seismic_hazard, 'disaggregation': psha_disagg}
        else:
            output['output']['psha'] = {"PGA": pga, "annual_rate_of_exceedance": seismic_hazard}
    if(output_plha):
        if(output_plha_disaggregation):
            for i in range(len(fsl)):
                plha_disagg[i] = plha_disagg[i] / liquefaction_hazard[i] * 100
            output['output']['plha'] = {"FSL": fsl, "annual_rate_of_nonexceedance": liquefaction_hazard, 'disaggregation': plha_disagg}
        else:
            output['output']['plha'] = {"FSL": fsl, "annual_rate_of_nonexceedance": liquefaction_hazard}
    return output
