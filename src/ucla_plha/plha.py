import numpy as np
import pandas as pd
import scipy as sp
from scipy.special import ndtr
import json
import liquefaction_models.cetin_et_al_2018 as cetin_et_al_2018
import liquefaction_models.idriss_boulanger_2012 as idriss_boulanger_2012
import ground_motion_models.ask14 as ask14
import ground_motion_models.bssa14 as bssa14
import ground_motion_models.ba14 as ba14
import ground_motion_models.cb14 as cb14
import geometry.geometry as geometry      

def get_source_data(source_type, source_model, dist_type, p_xyz, dist_cutoff):
    if(source_type == 'fault_source_models'):
        ruptures = pd.read_pickle('source_models/fault_source_models/' + source_model + '/ruptures.pkl', compression='gzip')
        fault_id = np.load('source_models/fault_source_models/' + source_model + '/fault_id.npy')
        m = ruptures['m'].values
        fault_type = ruptures['style'].values
        ruptures_segments = pd.read_pickle('source_models/fault_source_models/' + source_model + '/ruptures_segments.pkl', compression='gzip')
        segment_index = ruptures_segments['segment_index'].values
        if(dist_type == 'rjb'):
            tri = np.load('source_models/fault_source_models/' + source_model + '/tri_rjb.npy')
        elif(dist_type == 'rrup'):
            tri = np.load('source_models/fault_source_models/' + source_model + '/tri_rrup.npy')
        dist_all = geometry.point_triangle_distance(tri, p_xyz, fault_id)
        ruptures_segments['dist_all'] = dist_all[segment_index]
        dist = ruptures_segments.groupby('rupture_index')['dist_all'].min().values
        filter = (dist < dist_cutoff)
        rate = ruptures['rate'].values
    elif(source_type == 'point_source_models'):
        ruptures = pd.read_pickle('source_models/point_source_models/' + source_model + '/ruptures.pkl', compression='gzip')
        m = ruptures['m'].values
        fault_type = ruptures['style'].values
        node_index = np.load('source_models/point_source_models/' + source_model + '/node_index.npy')
        points = np.load('source_models/point_source_models/' + source_model + '/points.npy')
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
        mu_ln_pga, sigma_ln_pga = ba14.get_im(vs30, dist, m, fault_type)
    elif(gmm == 'ask14'):
        mu_ln_pga, sigma_ln_pga = ask14.get_im(vs30, dist, m, fault_type)
    else:
        print('incorrect ground motion model')
    return([mu_ln_pga, sigma_ln_pga])

def get_liquefaction_hazards(m, mu_ln_pga, sigma_ln_pga, rate, fsl, liquefaction_model, config):
    if(liquefaction_model=='cetin_et_al_2018'):
        c = config['liquefaction_models']['cetin_et_al_2018']
        fsl_hazards = cetin_et_al_2018.get_fsl_hazards(mu_ln_pga, sigma_ln_pga, rate, fsl, m, c['sigmav'], c['sigmavp'], c['vs12'], c['depth'], c['n160'], c['fc'], c['pa'])
    elif(liquefaction_model=='idriss_boulanger_2012'):
        c = config['liquefaction_models']['idriss_boulanger_2012']
        fsl_hazards = idriss_boulanger_2012.get_fsl_hazards(mu_ln_pga, sigma_ln_pga, rate, fsl, m, c['sigmav'], c['sigmavp'], c['depth'], c['n160'], c['fc'], c['pa'])
    return fsl_hazards

def get_hazard(config_file):
    # Validate config_file against schema
    #####
    # Write validation code here
    #####

    dist_types = {'bssa14': 'rjb', 'ba14': 'rjb', 'cb14': 'rjb', 'as14': 'rjb'}
    config = json.load(open(config_file))
    
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
            magnitude_bin_edges = np.asarray(config['output']['psha']['disaggregation']['magnitude_bin_edges'], dtype=float)
            distance_bin_edges = np.asarray(config['output']['psha']['disaggregation']['distance_bin_edges'], dtype=float)
            epsilon_bin_edges = np.asarray(config['output']['psha']['disaggregation']['epsilon_bin_edges'], dtype=float)
            magnitude_bin_center = 0.5 * (magnitude_bin_edges[0:-1] + magnitude_bin_edges[1:])
            distance_bin_center = 0.5 * (distance_bin_edges[0:-1] + distance_bin_edges[1:])
            epsilon_bin_center = 0.5 * (epsilon_bin_edges[0:-1] + epsilon_bin_edges[1:])
            psha_disagg = np.zeros((len(pga), len(magnitude_bin_center), len(distance_bin_center), len(epsilon_bin_center)))
        else:
            output_psha_disaggregation = False
        seismic_hazard = np.zeros(len(pga))
    else:
        output_psha = False
        output_psha_disaggregation = False

    if(config['output']['plha']):
        fsl = np.asarray(config['output']['plha']['fsl'], dtype=float)
        output_plha = True
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
                            for j in range(len(magnitude_bin_center)):
                                for k in range(len(distance_bin_center)):
                                    for l in range(len(epsilon_bin_center)):
                                        psha_disagg[i, j, k, l] += source_model_weight * ground_motion_model_weight * np.sum(seismic_hazards[i, (m >= magnitude_bin_edges[j]) & (m < magnitude_bin_edges[j+1]) & (rjb >= distance_bin_edges[k]) & (rjb < distance_bin_edges[k+1]) & (eps >= epsilon_bin_edges[l]) & (eps < epsilon_bin_edges[l+1])])
                        

                # Compute liquefaction hazard if requested in config file
                for liquefaction_model in config['liquefaction_models'].keys():
                    liquefaction_model_weight = config['liquefaction_models'][liquefaction_model]['weight']
                    if(output_plha):
                        liquefaction_hazards = get_liquefaction_hazards(m, mu_ln_pga, sigma_ln_pga, rate, fsl, liquefaction_model, config)
                        liquefaction_hazard += source_model_weight * ground_motion_model_weight * liquefaction_model_weight * np.sum(liquefaction_hazards, axis=0)

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
        output['output']["plha"] = {"FSL": fsl, "annual_rate_of_nonexceedance": liquefaction_hazard}
    return output