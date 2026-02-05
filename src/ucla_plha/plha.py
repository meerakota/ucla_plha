import json
import os
from importlib.resources import files

import jsonschema
import numpy as np
import pandas as pd
import scipy as sp
from scipy.stats import norm
from scipy.special import ndtr

from ucla_plha.liquefaction_models import (
    boulanger_idriss_2012,
    cetin_et_al_2018,
    ngl_smt_2024,
    moss_et_al_2006,
    boulanger_idriss_2016,
)
from ucla_plha.ground_motion_models import ask14, bssa14, cb14, cy14
from ucla_plha.geometry import geometry


def decompress_ucerf3_source_data():
    '''Decompresses ruptures.gz and ruptures_segments.gz in source_models/fault_source_models/ucerf3_fm31
    and source_models/fault_source_models/ucerf3_fm32. These files are compressed to reduce the
    size of the ucla_plha package, but decompressing the files each time the code is run is inefficient.
    This function decompresses them within the installation directory. The get_source_data() function
    checks for .pkl versions of the files, and uses them if they exist. Otherwise it uses the 
    .gz versions of the files and decompresses them at runtime. This function needs to be run only
    once.
    '''
    branches = ["ucerf3_fm31", "ucerf3_fm32"]
    for branch in branches:
        path = files("ucla_plha").joinpath(
            "source_models/fault_source_models/" + branch
        )
        ruptures = pd.read_pickle(str(path.joinpath("ruptures.gz")), compression="gzip")
        pd.to_pickle(ruptures, str(path.joinpath("ruptures.pkl")))
        ruptures_segments = pd.read_pickle(
            str(path.joinpath("ruptures_segments.gz")), compression="gzip"
        )
        pd.to_pickle(ruptures_segments, str(path.joinpath("ruptures_segments.pkl")))


def get_source_data(source_type, source_model, p_xyz, dist_cutoff, m_min, gmms):
    '''Returns magnitude, fault type, rate, distance, and fault geometry terms.
    
    Args:
        source_type (string): Either "fault_source_models" or "point_source_models"
        source_model (string): Directory for source_model within the source_type directory.
            Currently either "ucerf3_fm31" or "ucerf3_fm32"
        p_xyz (numpy array, dtype=float): Array containing x, y, z coordinates for point of interest, length = 3
        dist_cutoff (float): maximum distance to consider in seismic hazard analysis
        m_min (float): minimum magnitude to consider in seismic hazard analysis
        gmms (array, dtype=string): An array of strings defining ground motion models to use in seismic
            hazard analysis. Currently one or more of "ask14", "bssa14", "cb14", "cy14"

    Returns: A tuple containing the following arrays
        m (array, dtype=float): Numpy array of magnitudes, length = N
        fault_type (array, dtype=int): Numpy array of fault types. 1 = reverse, 2 = normal, 3 = strike slip, length = N
        rate (array, dtype=float): Numpy array of the rate of occurrence of each event, length = N
        rjb (array, dtype=float): Numpy array of Joyner-Boore distances between the site and each rupture, length = N
        rrup (array, dtype=float): Numpy array of rupture distance between the site and each rupture, length = N
        rx (array, dtype=float): Numpy array of distance from site to the surface projection of the top of each
            rupture, measured perpendicular to the strike, length = N
        rx1 (array, dtype=float): Numpy array of distance from site to the surface projection of the bottom of each
            rupture, measured perpendicular to the strike, length = N
        ry0 (array, dtype=float): Numpy array of the distance from the site to the surface projection of each
            rupture, measured parallel to the strike, length = N
        dip (array, dtype=float): Numpy array of dip angle for each rupture in degrees, length = N
        ztor (array, dtype=float): Numpy array of depth to the top of each rupture in km, length = N
        zbor (array, dtype=float): Numpy array of depth to the bottom of each rupture in km, length = N

    Notes:
        N = number of events
    '''
    if source_type == "fault_source_models":
        # Read files required by all ground motion models
        path = files("ucla_plha").joinpath(
            "source_models/fault_source_models/" + source_model
        )
        # Read decompressed version of files if they exist. Otherwise read zipped version.
        if os.path.exists(str(path.joinpath("ruptures.pkl"))):
            ruptures = pd.read_pickle(str(path.joinpath("ruptures.pkl")))
        else:
            ruptures = pd.read_pickle(
                str(path.joinpath("ruptures.gz")), compression="gzip"
            )
        m = ruptures["m"].values
        fault_type = ruptures["style"].values
        if os.path.exists(str(path.joinpath("ruptures_segments.pkl"))):
            ruptures_segments = pd.read_pickle(
                str(path.joinpath("ruptures_segments.pkl"))
            )
        else:
            ruptures_segments = pd.read_pickle(
                str(path.joinpath("ruptures_segments.gz")), compression="gzip"
            )
        segment_index = ruptures_segments["segment_index"].values
        rate = ruptures["rate"].values
        dip = ruptures["dip"].values
        ztor = ruptures["ztor"].values
        zbor = ruptures["zbor"].values

        # Now read files required by ask14, bssa14, cb14, and / or cy14
        # bssa14: rjb
        # ask14: rrup,rx,rx1,ry0,
        # cb14: rjb,rrup,rx
        # cy14: rjb,rrup,rx
        empty_array = np.empty(len(m))
        if any(gmm in ["ask14", "bssa14", "cb14", "cy14"] for gmm in gmms):
            tri_segment_id = np.load(str(path.joinpath("tri_segment_id.npy")))
        if any(gmm in ["bssa14", "cb14", "cy14"] for gmm in gmms):
            tri_rjb = np.load(str(path.joinpath("tri_rjb.npy")))
            rjb_all = geometry.point_triangle_distance(tri_rjb, p_xyz, tri_segment_id)
            ruptures_segments["rjb_all"] = rjb_all[segment_index]
        if any(gmm in ["ask14", "cb14", "cy14"] for gmm in gmms):
            rect_segment_id = np.load(str(path.joinpath("rect_segment_id.npy")))
            tri_rrup = np.load(str(path.joinpath("tri_rrup.npy")))
            rect = np.load(str(path.joinpath("rect_rjb.npy")))
            rrup_all = geometry.point_triangle_distance(tri_rrup, p_xyz, tri_segment_id)
            ruptures_segments["rrup_all"] = rrup_all[segment_index]
            Rx_all, Rx1_all, Ry0_all = geometry.get_Rx_Rx1_Ry0(
                rect, p_xyz, rect_segment_id
            )
            ruptures_segments["Rx_all"] = Rx_all[segment_index]
            ruptures_segments["Rx1_all"] = Rx1_all[segment_index]
            ruptures_segments["Ry0_all"] = Ry0_all[segment_index]
            ruptures_segments["dip_all"] = dip[segment_index]
            ruptures_segments["ztor_all"] = ztor[segment_index]
            ruptures_segments["zbor_all"] = zbor[segment_index]
        grouped_ruptures_segments = ruptures_segments.groupby("rupture_index")
        if any(gmm in ["ask14", "cb14", "cy14"] for gmm in gmms):
            rrup = grouped_ruptures_segments["rrup_all"].min().values
            filter = (rrup < dist_cutoff) & (m >= m_min)
            rx = grouped_ruptures_segments["Rx_all"].min().values
            rx1 = grouped_ruptures_segments["Rx1_all"].min().values
            ry0 = grouped_ruptures_segments["Ry0_all"].min().values
        else:
            rrup = empty_array
            rx = empty_array
            rx1 = empty_array
            ry0 = empty_array

        if any(gmm in ["bssa14", "cb14", "cy14"] for gmm in gmms):
            rjb = grouped_ruptures_segments["rjb_all"].min().values
            filter = (rjb < dist_cutoff) & (m >= m_min)
        else:
            rjb = empty_array

        return (
            m[filter],
            fault_type[filter],
            rate[filter],
            rjb[filter],
            rrup[filter],
            rx[filter],
            rx1[filter],
            ry0[filter],
            dip[filter],
            ztor[filter],
            zbor[filter],
        )

    elif source_type == "point_source_models":
        path = files("ucla_plha").joinpath(
            "source_models/point_source_models/" + source_model
        )
        ruptures = pd.read_pickle(
            str(path.joinpath("ruptures.pkl")), compression="gzip"
        )
        rate = ruptures["rate"].values
        m = ruptures["m"].values
        fault_type = ruptures["style"].values
        node_index = np.load(str(path.joinpath("node_index.npy")))
        points = np.load(str(path.joinpath("points.npy")))
        
        repi_temp = np.empty(np.max(node_index) + 1, dtype=float)
        for i, ni in enumerate(node_index):
            repi_temp[ni] = np.sqrt(
                (points[i, 0] - p_xyz[0]) ** 2
                + (points[i, 1] - p_xyz[1]) ** 2
                + (points[i, 2] - p_xyz[2]) ** 2
            )
        repi = repi_temp[ruptures["node_index"].values]
        rjb = repi / (1.0 + np.exp(-1.05 * (np.log(repi) - 1.037 * m + 4.2776)))

        filter = (rjb < dist_cutoff) & (m >= m_min)
        dip = np.empty(len(m), dtype=float)
        
        # using Kaklamanos et al. 2011 guidance for unknown dip, ztor, and zbor
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
        # Expression for Rrup is fit to default parameters (first row in table 1) in Thompson and Worden (2017)
        d = np.radians(80)
        zbor = ztor + w * np.sin(dip * np.pi / 180.0)
        rrup1 = np.sqrt((repi - 0.6 * w * np.cos(d)) ** 2 + ztor ** 2)
        filt = repi < 0.6 * w * np.cos(d) - ztor * np.tan(d)
        rrup1[filt] = (zhyp[filt] - repi[filt] * np.tan(d)) / (np.cos(d) + np.sin(d) * np.tan(d))
        rrup2 = zhyp
        filt = repi > 1.7 * w
        rrup2[filt] = np.sqrt((repi[filt] - 1.7 * w[filt] / 2.0) ** 2 + zhyp[filt] ** 2)
        rrup = 0.5 * rrup1 + 0.5 * rrup2
        # Some assumptions here for rx. The rrup model is dominated by the site being updip from fault, so make rx positive.
        rx = -rjb / np.sqrt(2.0)
        rx[rjb == 0] = 0.5 * w[rjb == 0] * np.cos(d)
        # Calculate rx1 after computing rx
        rx1 = rx + w * np.cos(d)
        # Compute ry0 the same way as rx, but use fault length instead of width. Assume aspect ratio of 1.7
        ry0 = rjb / np.sqrt(2.0)
        ry0[rjb==0] = 0.5 * 1.7 * w[rjb == 0] * np.cos(d)
        
        return (
            m[filter],
            fault_type[filter],
            rate[filter],
            rjb[filter],
            rrup[filter],
            rx[filter],
            rx1[filter],
            ry0[filter],
            dip[filter],
            ztor[filter],
            zbor[filter],
        )


def get_ground_motion_data(
    gmm,
    vs30,
    measured_vs30,
    z1p0,
    z2p5,
    fault_type,
    rjb,
    rrup,
    rx,
    rx1,
    ry0,
    m,
    ztor,
    zbor,
    dip,
):
    '''Computes arrays containing mean and standard deviation of the natural log of a ground
    motion intensity measure.

    Inputs:
        gmm (string): Ground motion model. One of "ask14", "bssa14", "cb14", "cy14"
        vs30 (float): Time-averaged shear wave velocity in the upper 30m in m/s
        measured_vs30 (bool): boolean field indicating whether vs30 is measured (True) or inferred (False)
        z1p0 (float): Isosurface depth to a shear wave velocity of 1.0 km/s in km, length = N
        z2p5 (float): Isosurface depth to a shear wave velocity of 2.5 km/s in km, length = N
        fault_type (Numpy array, dtype=int): Numpy array of fault type. 1 = reverse, 2 = normal, 3 = strike slip, length = N
        rjb (Numpy array, dtype=float): Numpy array of Joyner-Boore distances between the site and each rupture, length = N
        rrup (Numpy array, dtype=float): Numpy array of rupture distance between the site and each rupture, length = N
        rx (Numpy array, dtype=float): Numpy array of distance from site to the surface projection of the top of each 
            rupture, measured perpendicular to the strike, length = N
        rx1 (Numpy array, dtype=float): Numpy array of distance from site to the surface projection of the bottom of each
            rupture, measured perpendicular to the strike, length = N
        ry0 (Numpy array, dtype=float): Numpy array of the distance from the site to the surface projection of each
            rupture, measured parallel to the strike, length = N
        m (Numpy array, dtype=float): Numpy array of magnitudes, length = N
        ztor (Numpy array, dtype=float): Numpy array of depth to the top of each rupture in km, length = N
        zbor (Numpy array, dtype=float): Numpy array of depth to the bottom of each rupture in km, length = N
        dip (Numpy array, dtype=float): Numpy array of dip angle for each rupture in degrees, length = N

    Returns:
        mu_ln_pga (Numpy array, dtype=float): array of the mean of the natural logs of the ground motion intensity measure, IM
        sigma_ln_pga (Numpy array, dtype=float): array of the standard deviation of the natural logs of the IM
    
    Note:
        N = number of events
    '''
    if gmm == "bssa14":
        mu_ln_pga, sigma_ln_pga = bssa14.get_im(vs30, rjb, m, fault_type)
    elif gmm == "cb14":
        mu_ln_pga, sigma_ln_pga = cb14.get_im(
            vs30, rjb, rrup, rx, rx1, m, fault_type, ztor, zbor, dip, z2p5=z2p5
        )
    elif gmm == "cy14":
        mu_ln_pga, sigma_ln_pga = cy14.get_im(
            vs30, rjb, rrup, rx, m, fault_type, measured_vs30, dip, ztor, z1p0=z1p0
        )
    elif gmm == "ask14":
        mu_ln_pga, sigma_ln_pga = ask14.get_im(
            vs30, rrup, rx, rx1, ry0, m, fault_type, measured_vs30, dip, ztor, z1p0=z1p0
        )
    else:
        print("incorrect ground motion model")
    return [mu_ln_pga, sigma_ln_pga]


def get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config):
    '''Computes log-normal cumulative distribution functions for factor of safety of liquefaction

    Inputs:
        m (array, dtype=float): Numpy array of magnitudes, length = N
        mu_ln_pga (array, dtype=float): Numpy array of the mean of the natural logs of the earthquake
            ground motion intensity measure, length = N
        sigma_ln_pga (array, dtype=float): Numpy array of the standard deviation of the natural logs
            of the ground motion intensity measure, length = N
        fsl (array, dtype=float): Numpy array of factor of safety values at which to compute the 
            liquefaction hazard curve, length = L
        liquefaction_model (string): Liquefaction model to use. One of "cetin_et_al_2018", "moss_et_al_2006",
            "boulanger_idriss_2016", "boulanger_idriss_2012", "ngl_smt_2024".
        config (dict): A Python dictionary read from the config file
    
    Returns:
        fsl_cdfs (Numpy ndarray, dtype=float): Numpy array of cumulative distribution functions representing
            down-crossing rate of fsl for each event, length = N x L
        eps (Numpy ndarray, dtype=float): Numpy array of epsilon values, representing number of standard deviations
            of the natural log of fsl relative to the mean of the natural log of fsl, length = N x L
        
    Notes:
        N = number of events
        L = number of fsl values at which to compute hazard
    '''
    if liquefaction_model == "cetin_et_al_2018":
        c = config["liquefaction_models"]["cetin_et_al_2018"]
        return cetin_et_al_2018.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["vs12"],
            c["depth"],
            c["n160"],
            c["fc"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "moss_et_al_2006":
        c = config["liquefaction_models"]["moss_et_al_2006"]
        return moss_et_al_2006.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["qc"],
            c["fs"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "boulanger_idriss_2016":
        c = config["liquefaction_models"]["boulanger_idriss_2016"]
        return boulanger_idriss_2016.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["qc1ncs"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "boulanger_idriss_2012":
        c = config["liquefaction_models"]["boulanger_idriss_2012"]
        return boulanger_idriss_2012.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            c["sigmav"],
            c["sigmavp"],
            c["depth"],
            c["n160"],
            c["fc"],
            fsl,
            c["pa"],
        )
    elif liquefaction_model == "ngl_smt_2024":
        c = config["liquefaction_models"]["ngl_smt_2024"]
        return ngl_smt_2024.get_fsl_cdfs(
            mu_ln_pga,
            sigma_ln_pga,
            m,
            np.asarray(c["ztop"], dtype=float),
            np.asarray(c["zbot"], dtype=float),
            np.asarray(c["qc1ncs"], dtype=float),
            np.asarray(c["ic"], dtype=float),
            np.asarray(c["sigmav"], dtype=float),
            np.asarray(c["sigmavp"], dtype=float),
            np.asarray(c["ksat"], dtype=float),
            fsl,
            float(c.get("pa", 101.325)),
        )


def get_disagg(hazards, m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges):
    """Computes disaggregation for seismic hazard and/or liquefaction hazard

    Inputs:
        hazards (Numpy ndarray, dtype = float) = array of hazard values, shape = L x N.
        m (Numpy array, dtype = float) = array of magnitude values, length = N
        r (Numpy array, dtype = float) = array of distance values, length = N
        eps (Numpy ndarray, dtype = float) = array of epsilon values, shape = L x N
        m_bin_edges (Numpy array, dtype = float) = array defining magnitude bin edges, length = M + 1
        r_bin_edges (Numpy array, dtype = float) = array defining distance bin edges, length = R + 1
        eps_bin_edges (Numpy array, dtype = float) = array defining epsilon bin edges, length = E + 1
    
    Returns:
        disagg (Numpy ndarray, dtype=float) = Numpy array of contribution to hazard within each 
            magnitude, distance, and epsilon bin for each pga (PSHA) or fsl (PLHA) value, shape = L x M x R x E

    Notes:
        N = number of events
        L = number of intensity measure values
        M = number of magnitude bins (note that m_bin_edges has a length of M + 1)
        R = number of distance bins (note that r_bin_edges has a length of R + 1)
        E = number of epsilon bins (note that eps_bin_edges has a length of E + 1)
    """
    # use Numpy digitize function to assign bin numbers
    m_hazard = np.digitize(m, m_bin_edges)
    r_hazard = np.digitize(r, r_bin_edges)
    eps_hazard = np.digitize(eps, eps_bin_edges)
    
    # compute number of bins per intensity measure value
    Nbins = (len(m_bin_edges) - 1) * (len(r_bin_edges) - 1) * (len(eps_bin_edges) - 1)

    # use tensor rank reduction to create L x Nbins length array of indices
    bin_indices = (
        eps_hazard
        - 1
        + (r_hazard - 1) * (len(eps_bin_edges) - 1)
        + (m_hazard - 1) * (len(eps_bin_edges) - 1) * (len(r_bin_edges) - 1)
    )
    
    # create empty array to store hazard sums, and use Numpy bincount to efficiently sum hazards
    # within each bin. The np.bincount function must be performed on a 1D array, so we need to loop
    # over the pga values. The cost of that loop is minimal since the number of pga values
    # is generally small
    sum_hazard = np.empty((len(bin_indices), Nbins))
    for i in range(len(bin_indices)):
        sum_hazard[i] = np.bincount(bin_indices[i], weights=hazards[i], minlength=Nbins)

    # reshape the sum_hazard array
    disagg = sum_hazard.reshape(
        (len(bin_indices), len(m_bin_edges) -1, len(r_bin_edges) - 1, len(eps_bin_edges) - 1)
    )

    return disagg


def get_hazard(config_file):
    '''Reads config file and runs PSHA and PLHA

    Inputs:
        config_file (string): Filename, including path, of config file. Must follow the schema
            defined by ucla_plha_schema.json. See documentation for more thorough documentation
            of the config file.
    
    Returns:
        output (dict): Python dictionary containing output of analysis. The output contains all
            of the inputs for preservation, along with the hazard curve(s) and any requested
            disaggregation data. See documentation for more thorough description of output.
    '''
    # Validate config_file against schema. If ngl_smt_2024 liquefaction model is used, the cpt_data file is
    # validated in the get_liquefaction_hazards function.

    schema = json.loads(open(files("ucla_plha").joinpath("ucla_plha_schema.json")).read().lower())
    config = json.loads(open(config_file).read().lower())

    # validate config file and return messageg if errors are encountered
    try:
        jsonschema.validate(config, schema)
    except jsonschema.ValidationError as e:
        print("Config File Error:", e.message)
        return
    
    # normalize weights in config file
    fault_source_model_weight_sum = 0.0
    fault_source_models = ['ucerf3_fm31', 'ucerf3_fm32']
    for fault_source_model in fault_source_models:
        fault_source_model_weight_sum += config['source_models'].get('fault_source_models', {}).get(fault_source_model, {}).get('weight', 0.0)
    if(fault_source_model_weight_sum > 0):
        for fault_source_model in fault_source_models:
            if(config['source_models'].get('fault_source_models', {}).get(fault_source_model, {})):
                config['source_models'].get('fault_source_models', {})[fault_source_model]['weight'] /= fault_source_model_weight_sum

    point_source_model_weight_sum = 0.0
    point_source_models = ['ucerf3_fm31_grid_sub_seis', 'ucerf3_fm31_grid_unassociated', 'ucerf3_fm32_grid_sub_seis', 'ucerf3_fm32_grid_unassociated']
    for point_source_model in point_source_models:
        point_source_model_weight_sum += config['source_models'].get('point_source_models', {}).get(point_source_model, {}).get('weight', 0.0)
    if(point_source_model_weight_sum > 0):
        for point_source_model in point_source_models:
            if(config['source_models']['point_source_models'].get(point_source_model, {})):
                config['source_models']['point_source_models'][point_source_model]['weight'] /= point_source_model_weight_sum

    ground_motion_model_weight_sum = 0.0
    ground_motion_models = ['bssa14', 'ask14', 'cb14', 'cy14']
    for ground_motion_model in ground_motion_models:
        ground_motion_model_weight_sum += config['ground_motion_models'].get(ground_motion_model, {}).get('weight', 0.0)
    if(ground_motion_model_weight_sum > 0):
        for ground_motion_model in ground_motion_models:
            if(config['ground_motion_models'].get(ground_motion_model, {})):
                config['ground_motion_models'][ground_motion_model]['weight'] /= ground_motion_model_weight_sum

    liquefaction_model_weight_sum = 0.0
    liquefaction_models = ['boulanger_idriss_2012', 'boulanger_idriss_2016', 'cetin_et_al_2018', 'moss_et_al_2006', 'ngl_smt_2024']
    for liquefaction_model in liquefaction_models:
        liquefaction_model_weight_sum += config['liquefaction_models'].get(liquefaction_model, {}).get('weight', 0.0)
    if(liquefaction_model_weight_sum > 0):
        for liquefaction_model in liquefaction_models:
            if(config['liquefaction_models'].get(liquefaction_model, {})):
                config['liquefaction_models'][liquefaction_model]['weight'] /= liquefaction_model_weight_sum


    # Read site properties
    latitude = config["site"]["latitude"]
    longitude = config["site"]["longitude"]
    elevation = config["site"]["elevation"]
    point = np.asarray([latitude, longitude, elevation])
    dist_cutoff = config["site"]["dist_cutoff"]
    m_min = config["site"]["m_min"]
    p_xyz = geometry.point_to_xyz(point)

    # Read output properties
    if "psha" in config["output"].keys():
        pga = np.asarray(config["output"]["psha"]["pga"], dtype=float)
        output_psha = True
        if "disaggregation" in config["output"]["psha"].keys():
            output_psha_disaggregation = True
            psha_magnitude_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["magnitude_bin_edges"],
                dtype=float,
            )
            psha_distance_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["distance_bin_edges"],
                dtype=float,
            )
            psha_epsilon_bin_edges = np.asarray(
                config["output"]["psha"]["disaggregation"]["epsilon_bin_edges"],
                dtype=float,
            )
            psha_magnitude_bin_center = 0.5 * (
                psha_magnitude_bin_edges[0:-1] + psha_magnitude_bin_edges[1:]
            )
            psha_distance_bin_center = 0.5 * (
                psha_distance_bin_edges[0:-1] + psha_distance_bin_edges[1:]
            )
            psha_epsilon_bin_center = 0.5 * (
                psha_epsilon_bin_edges[0:-1] + psha_epsilon_bin_edges[1:]
            )
            psha_disagg = np.zeros(
                (
                    len(pga),
                    len(psha_magnitude_bin_center),
                    len(psha_distance_bin_center),
                    len(psha_epsilon_bin_center),
                ),
                dtype=float,
            )
        else:
            output_psha_disaggregation = False
        seismic_hazard = np.zeros(len(pga))
    else:
        output_psha = False
        output_psha_disaggregation = False

    if "plha" in config["output"].keys():
        fsl = np.asarray(config["output"]["plha"]["fsl"], dtype=float)
        output_plha = True
        if "disaggregation" in config["output"]["plha"].keys():
            output_plha_disaggregation = True
            plha_magnitude_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["magnitude_bin_edges"],
                dtype=float,
            )
            plha_distance_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["distance_bin_edges"],
                dtype=float,
            )
            plha_epsilon_bin_edges = np.asarray(
                config["output"]["plha"]["disaggregation"]["epsilon_bin_edges"],
                dtype=float,
            )
            plha_magnitude_bin_center = 0.5 * (
                plha_magnitude_bin_edges[0:-1] + plha_magnitude_bin_edges[1:]
            )
            plha_distance_bin_center = 0.5 * (
                plha_distance_bin_edges[0:-1] + plha_distance_bin_edges[1:]
            )
            plha_epsilon_bin_center = 0.5 * (
                plha_epsilon_bin_edges[0:-1] + plha_epsilon_bin_edges[1:]
            )
            plha_disagg = np.zeros(
                (
                    len(pga),
                    len(plha_magnitude_bin_center),
                    len(plha_distance_bin_center),
                    len(plha_epsilon_bin_center),
                ),
                dtype=float,
            )
        else:
            output_plha_disaggregation = False
        liquefaction_hazard = np.zeros(len(fsl))
    else:
        output_plha = False
        output_plha_disaggregation = False

    # Loop over all ground motion models to get list of distance types
    gmms = []
    for gmm in config["ground_motion_models"].keys():
        gmms.append(gmm)
    # Loop over source models. We have fault_source_models and point_source_models, so there are two loops
    for source_model in config["source_models"].keys():
        for fault_source_model in config["source_models"][source_model].keys():
            m, fault_type, rate, rjb, rrup, rx, rx1, ry0, dip, ztor, zbor = (
                get_source_data(
                    source_model, fault_source_model, p_xyz, dist_cutoff, m_min, gmms
                )
            )
            source_model_weight = config["source_models"][source_model][
                fault_source_model
            ]["weight"]
            # Loop over ground motion models.
            for ground_motion_model in config["ground_motion_models"].keys():
                # retrieve parameters common to all models
                ground_motion_model_weight = config["ground_motion_models"][
                    ground_motion_model
                ]["weight"]
                # move on to next ground motion model if weight is less than or equal to zero
                if (ground_motion_model_weight <= 0):
                    continue
                vs30 = config["ground_motion_models"][ground_motion_model]["vs30"]
                z1p0 = None
                z2p5 = None
                measured_vs30 = False
                # retrieve model-specific parameters
                if (ground_motion_model == "ask14") or (ground_motion_model == "cy14"):
                    if (
                        "measured_vs30"
                        in config["ground_motion_models"][ground_motion_model].keys()
                    ):
                        measured_vs30 = config["ground_motion_models"][ground_motion_model][
                            "measured_vs30"
                        ]
                    if (
                        "z1p0"
                        in config["ground_motion_models"][ground_motion_model].keys()
                    ):
                        z1p0 = config["ground_motion_models"][ground_motion_model][
                            "z1p0"
                        ]
                if ground_motion_model == "cb14":
                    if (
                        "z2p5"
                        in config["ground_motion_models"][ground_motion_model].keys()
                    ):
                        z2p5 = config["ground_motion_models"][ground_motion_model][
                            "z2p5"
                        ]
                mu_ln_pga, sigma_ln_pga = get_ground_motion_data(
                    ground_motion_model,
                    vs30,
                    measured_vs30,
                    z1p0,
                    z2p5,
                    fault_type,
                    rjb,
                    rrup,
                    rx,
                    rx1,
                    ry0,
                    m,
                    ztor,
                    zbor,
                    dip,
                )
                # Compute seismic hazard if requested in config file
                if output_psha:
                    eps = (np.log(pga[:, np.newaxis]) - mu_ln_pga) / sigma_ln_pga
                    seismic_hazards = (1 - ndtr(eps)) * rate
                    seismic_hazard += (
                        source_model_weight
                        * ground_motion_model_weight
                        * np.sum(seismic_hazards, axis=1)
                    )
                    # Compute seismic hazard disaggregation if requested in config file
                    if output_psha_disaggregation:
                        psha_disagg += (
                            source_model_weight
                            * ground_motion_model_weight
                            * get_disagg(
                                seismic_hazards,
                                m,
                                rjb,
                                eps,
                                psha_magnitude_bin_edges,
                                psha_distance_bin_edges,
                                psha_epsilon_bin_edges,
                            )
                        )
                # Compute liquefaction hazard if requested in config file
                if "liquefaction_models" in config.keys():
                    for liquefaction_model in config["liquefaction_models"].keys():
                        liquefaction_model_weight = config["liquefaction_models"][
                            liquefaction_model
                        ]["weight"]
                        if (liquefaction_model_weight <= 0):
                            continue
                        if output_plha:
                            liquefaction_hazards, eps = get_liquefaction_cdfs(
                                m,
                                mu_ln_pga,
                                sigma_ln_pga,
                                fsl,
                                liquefaction_model,
                                config,
                            )
                            liquefaction_hazards *= rate[:, np.newaxis]
                            liquefaction_hazard += (
                                source_model_weight
                                * ground_motion_model_weight
                                * liquefaction_model_weight
                                * np.sum(liquefaction_hazards, axis=0)
                            )
                            eps = eps.T
                            liquefaction_hazards = liquefaction_hazards.T
                            # Compute liquefaction hazard disaggregation if requested in config file
                            if output_plha_disaggregation:
                                plha_disagg += (
                                    source_model_weight
                                    * ground_motion_model_weight
                                    * liquefaction_model_weight
                                    * get_disagg(
                                        liquefaction_hazards,
                                        m,
                                        rjb,
                                        eps,
                                        plha_magnitude_bin_edges,
                                        plha_distance_bin_edges,
                                        plha_epsilon_bin_edges,
                                    )
                                )
    # Now prepare output
    output = {}
    output["input"] = config
    output["output"] = {}
    if output_psha:
        if output_psha_disaggregation:
            for i in range(len(pga)):
                psha_disagg[i] = psha_disagg[i] / seismic_hazard[i] * 100.0
            output["output"]["psha"] = {
                "PGA": pga.tolist(),
                "annual_rate_of_exceedance": seismic_hazard.tolist(),
                "disaggregation": psha_disagg.tolist(),
            }
        else:
            output["output"]["psha"] = {
                "PGA": pga.tolist(),
                "annual_rate_of_exceedance": seismic_hazard.tolist(),
            }
    if output_plha:
        if output_plha_disaggregation:
            for i in range(len(fsl)):
                plha_disagg[i] = plha_disagg[i] / liquefaction_hazard[i] * 100.0
            output["output"]["plha"] = {
                "FSL": fsl.tolist(),
                "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
                "disaggregation": plha_disagg.tolist(),
            }
        else:
            output["output"]["plha"] = {
                "FSL": fsl.tolist(),
                "annual_rate_of_nonexceedance": liquefaction_hazard.tolist(),
            }

    if "outputfile" in config["output"].keys():
        if config["output"]["outputfile"] == "default":
            outputfilename = config_file.split(".json")[0] + "_output.json"
        else:
            outputfilename = config["output"]["outputfile"]
        with open(outputfilename, "w") as outputfile:
            json.dump(output, outputfile, indent=4)

    return output


