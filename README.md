# ucla_plha: A Probabilistic Liquefaction Hazard Assessment Python Package

**Meera L. Kota <meerakota@g.ucla.edu> and Scott J. Brandenberg <sjbrandenberg@g.ucla.edu>**  

## Installation
> pip install ucla_plha

## Overview 
The ucla_plha Python package performs probabilistic liquefaction hazard analysis (PLHA) to quantify the annual down-crossing rate (or alternatively, the annual rate of non-exceedance) of factor of safety of liquefaction. The code integrates probabilistic seismic hazard analysis (PSHA) within the same computational engine as PLHA. The ucla_plha_code performs PLHA inside the hazard integral, which distinguishes it from other codes that decouple the PSHA from the PLHA. The benefits of performing PLHA inside the hazard integral are (1) it is computationally efficient when ground motion demand and liquefaction capacity are represented by probability density functions of the same functional form (i.e., log-normal) because convolution can be performed quickly in closed-form, (2) it inherently accounts for the influence of magnitude on liquefaction because liquefaction assessment is performed for every event, and (3) disaggregation can be performed on the liquefaction hazard analysis. The cons are (1) computer code must be developed to integrate PSHA with PLHA (which is the purpose of the ucla_plha package), and (2) convolution becomes computationally expensive when ground motion demand and liquefaction capacity for a given event have different distribution functions. A benefit of decoupling PSHA and PLHA is that the PSHA can be performed using any available seismic hazard code, as long as a disaggregation is provided. The downside is that users must interpret disaggregation data at multiple return periods to identify a weighted set of magnitudes for performing liquefaction hazard analysis. This step is time consuming, involves judgment, and is a numerical approximation that asymptotically approaches the exact solution as the number of considered magnitudes approaches the number of distinct events considered in the PSHA.  

Currently, the code uses the Uniform California Earthquake Rupture Forecast, Version 3 (UCERF3, Field et al. 2013) source model, and is therefore applicable for shallow crustal earthquakes in California. Source-to-site distances are computed by dividing fault segments into triangles, and using a vectorized version of the closed-form point-to-triangle distance algorithm developed by [Eberly 1999](https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf). Background seismicity is modeled as point sources in UCERF3, and we treat each point as the epicenter and compute $R_{EPI}$ and subsequently relate it to $R_{RUP}$ and $R_{JB}$ using an approximation of algorithm by [Thompson and Worden 2017](https://doi.org/10.1785/0120170174). The code uses the four NGAWest2 ground motion models that model site response using $V_{S30}$, including Abrahamson et al. (2014), Boore et al. (2014), Campbell and Bozorgnia (2014), and Chiou and Youngs (2014). The code currently uses five probabilistic liquefaction models, including two standard penetration test models (Boulanger and Idriss 2012, and Cetin et al. 2018) and three cone penetration test models (Moss et al. 2006, Boulanger and Idriss 2016, and Ulmer et al. 2024).

## Example Usage
> import ucla_plha  
> output = ucla_plha.get_hazard("config.json")

## Config File Format

The ucla_plha code uses a [Javascript Object Notation (JSON)](https://www.json.org/json-en.html) input file to configure the analysis, and uses [JSON Schema](https://json-schema.org/) to validate user inputs. Example config files are provided in the examples directory, and the schema is defined by ucla_plha_schema.json located in the src/ucla_plha directory. The config file is case insensitive..

The structure of the config file is provided below. Note that the keys are all lowercase in the example below, which follows the PEP 8 convention for Python variable names. It may seem counter-intuitive to use "n160" instead of "N160" because the PEP 8 convention does not align with geotechnical convention. We have opted to stick with PEP 8 convention here.

```json
{
    "site": {
        "latitude": latitude,
        "longitude": longitude,
        "elevation": elevation,
        "dist_cutoff": dist_cutoff,
        "m_min": m_min
    },
    "source_models": {
        "fault_source_models": {
            "ucerf3_fm31": {
                "weight": weight
            },
            "ucerf3_fm32": {
                "weight": weight
            }
        },
        "point_source_models": {
            "ucerf3_fm31": {
                "weight": weight
            },
            "ucerf3_fm32": {
                "weight": weight
            }
        }
    },
    "ground_motion_models": {
        "ask14": {
            "vs30": vs30,
            "measured_vs30": measured_vs30,
            "z1p0": z1p0,
            "weight": weight
        },
        "bssa14": {
            "vs30": vs30,
            "weight": weight
        },
        "cb14": {
            "vs30": vs30,
            "measured_vs30": measured_vs30,
            "z2p5": z2p5,
            "weight": weight
        },
        "cy14": {
            "vs30": vs30,
            "measured_vs30": measured_vs30,
            "weight": weight
        }
    },
    "liquefaction_models": {
        "boulanger_idriss_2012": {
            "n160": n160,
            "fc": fc,
            "sigmav": sigmav,
            "sigmavp": sigmavp,
            "depth": depth,
            "pa": pa,
            "weight": weight
        },
        "boulanger_idriss_2016": {
            "qc1Ncs": qc1Ncs,
            "sigmav": sigmav,
            "sigmavp": sigmavp,
            "depth": depth,
            "pa": pa,
            "weight": weight
        },
        "cetin_et_al_2018": {
            "n160": n160,
            "fc": fc,
            "vs12": vs12,
            "sigmav": sigmav,
            "sigmavp": sigmavp,
            "depth": depth,
            "pa": pa,
            "weight": weight
        },
        "moss_et_al_2006": {
            "qc": qc,
            "fs": fs,
            "sigmav": sigmav,
            "sigmavp": sigmavp,
            "depth": depth,
            "pa": pa,
            "weight": weight
        },
        "ngl_smt_2024": {
            "ztop": [ztop_array],
            "zbot": [zbot_array],
            "qc1Ncs": [qc1Ncs_array],
            "Ic": [Ic_array],
            "sigmav": [sigmav_array],
            "sigmavp": [sigmavp_array],
            "Ksat": [Ksat_array],
            "pa": pa,
            "weight": weight
        }
    },
    "output": {
        "psha": {
            "pga": [pga_array],
            "disaggregation": {
                "magnitude_bin_edges": [magnitude_bin_edges_array],
                "distance_bin_edges": [distance_bin_edges_array],
                "epsilon_bin_edges": [espilon_bin_edges_array]
            }
        },
        "plha": {
            "fsl": [fsl_array],
            "disaggregation": {
                "magnitude_bin_edges": [magnitude_bin_edges_array],
                "distance_bin_edges": [distance_bin_edges_array],
                "epsilon_bin_edges": [espilon_bin_edges_array]
            }
        },
        "outputfile": outputfile
    }
}
```

### Top Level Keys

The top level keys are all JSON objects. If the "liquefaction_models" key and its associated JSON objects are excluded from the config file, the code will perform only a PSHA without performing a PLHA.  

| key | required |
| --- | ---------- |
| "site" | x|
| "source_models" | x |
| "ground_motion_models" | x |
| "liquefaction_models" |  |
| "output" | x |

### "site" Keys  

| key | definition | units | required |
| --- | ---------- | ----- | -------- |
| "latitude" | latitude | decimal degrees | x |
| "longitude" | longitude | decimal degrees | x |
| "elevation" | elevation | km | x |
| "dist_cutoff" | maximum distance to use in PSHA | km |   |
| "m_min" | minimum magnitude to use in PSHA | - | |

### "source_models" Keys

The "source_model" JSON object is not well suited to representation in tabular form because it has a high nesting level, but only a single variable. The 2nd nesting level contains two keys "fault_source_models" and "point_source_models". The 3rd nesting level contains two keys "ucerf3_fm31" and "ucerf3_fm32". The 4th nesting level contains a single key "weight". The weights are applied to the hazard curves computed using each source model, and summed together. The weights within the "fault_source_models" object should sum to unity, and the weights within the "point_source_models" object should sum to unity. The weights can be input however the user sees fit into the config file. However, the weights will then be normalized (divided by the sum total) within the code. This holds true for the fault source models, point source models, ground motion models, and liquefaction models. 

### "ground_motion_models" Keys

The 2nd nesting level keys define the ground motion model, and must be one of "ask14", "bssa14", "cb14", and/or "cy14". The 3rd nesting level keys define the variables indicated in the table below. 

| key | definition | units | notes |
| --- | ---------- | ----- | ----- |
| "vs30" | time-averaged shear wave velocity in upper 30m. | m/s | required for all models |
| "measured_vs30" | boolean indicator of whether vs30 was measured or inferred | - | optional for "ask14" and "cy14". default = false |
| "z1p0" | depth to isosurface where $V_S$ = 1.0 km/s | km | optional for "ask14" and "cy14" |
| "z2p5" | depth to isosurface where $V_S$ = 2.5 km/s | km | optional for "cb14" |
| "weight" | weight to assign to ground motion model | - | required for all models |

Notes: 
* $z_{1.0}$ is used in the BSSA14 model, but has no influence on PGA. For that reason, it is not allowed as an input here. 
* Isosurface depths $z_{1.0}$ and $z_{2.5}$ are optional. When they are not specified, the "centered" value corresponding to the specified $V_{S30}$ is used. This is equivalent to using $\delta z$ = 0.

### "liquefaction_model" Keys

The 2nd nesting level keys define the liquefaction model, and must be one of "boulanger_idriss_2012", "boulanger_idriss_2016", "cetin_et_al_2018", "moss_et_al_2006", "ngl_smt_2024". The 3rd nesting level keys define the variables used by the liquefaction model.

| key | definition | units | notes |
| --- | ---------- | ----- | ----- |
| "n160" | overburden- and energy-corrected SPT blow count | - | required for "boulanger_idriss_2012" and "cetin_et_al_2018" |
| "fc" | fines content | % | required for "boulanger_idriss_2012" and "cetin_et_al_2018" |
| "vs12" | time-averaged shear wave velocity in upper 12 m | m/s | required for "cetin_et_al_2018" |
| "qc1Ncs" | overburden- and fines-corrected dimensionless cone tip resistance | - | required for "boulanger_and_idriss_2016" and "ngl_smt_2024" |
| "qc" | cone tip resistance | MPa | required for "moss_et_al_2006" |
| "fs" | cone sleeve friction | MPa | required for "moss_et_al_2006" |
| "Ic" | soil behavior type index | - | required for "ngl_smt_2024" |
| "Ksat" | saturation factor applied to probability of triggering | - | required for "ngl_smt_2024" |
| "sigmav" | vertical total stress at center of layer | kPa | required for all models |
| "sigmavp" | vertical effective stress at center of layer | kPa | required for all models |
| "depth" | depth to center of layer | m | required for all models |
| "pa" | atmospheric pressure | kPa | required for all models |
| "weight" | weight to assign to liquefaction model | - | required for all models |

Note, you cannot include keys if they are not required by a liquefaction model. For example, you cannot include "n160" in one of the cone penetration test models.

### "output" Keys

The "output" JSON object specifies whether to output PSHA and / or PLHA results, and if so, whether to output disaggregation data. Disaggregation is computationally demanding. If you do not want disaggregation data, simply leave out the "disaggregation" key and associated objects from your config file, and the analysis will run more quickly.

| key | definition | notes |
| --- | ---------- | ----- |
| "psha" | include this key to output a PSHA hazard curve | optional |
| "plha" | include this key to output a PLHA hazard curve | optional |
| "disaggregation" | include this key under "psha" or "plha" to output disaggregation data | optional |
| "magnitude_bin_edges" | edges of magnitude bins for disaggregation | required if "disaggregation" is included |
| "distance_bin_edges" | edges of distance bins for disaggregation | required if "disaggregation" is included |
| "epsilon_bin_edges" |edges of epsilon bins for disaggregation | required if "disaggregation" is included |
| "outputfile" | name of file for saving output | optional |

Notes:
* You can specify the "plha" key and its "disaggregation" object without also including the "psha" key. In this case, the PSHA will still be performed, but the results will not be output. 
* Whether or not you specify the "outputfile" key, the code will return a Python dictionary containing the data.
* The input JSON object will be saved to the outputfile so that the inputs and outputs are always together. This will prevent inconsistent results if, for example, an analysis is run on an input file that is subsequently modified. Simply saving the input and output files in the same directory does not guard against this inconsistency.


## Functions

Functions available in the ucla_plha package are divided into subpackages
| function | description |
| -------- | ----------- |
| [```decompress_ucerf3_source_data()```](#function1kwargs) | Parse and compress UCERF3 source code   |
| [```get_source_data(source_type, source_model, p_xyz, dist_cutoff, gmms)```](#function2args) | Read files required for Ground Motion Models including distances |
| [```get_ground_motion_data(gmm, vs30, fault_type, rjb, rrup, rx, rx1, ry0, m, ztor, zbor, dip, z1p0, z2p5, measured_vs30)```](#function3args) | Function to run the 4 $V_{S30}$ dependent ground motion models|
| [```get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config)```](#function3args) | Running probabilistic liquefaction triggering models |
| [```get_disagg(hazards, m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges)```](#function3args) | Running seismic and liquefaction disaggregation |
| [```get_hazard(config_file)```](#function3args) | Run hazard calculation using ground motion and probabilistic triggering model outputs |

## decompress_ucerf3_source_data() 

Some of the input files for the UCERF3 model are quite large and are compressed to facilitate more efficient package installation. Specifically, within ucerf3_fm31 and ucerf3_fm32 directories in the source_models/fault_source_models, the ruptures.gz, and ruptures_segments.gz are zipped pickle (.pkl) files. Using the compressed files in the ucla_plha code results is inefficient becuase the files must be decompressed each time the code is run. This function decompresses the ruptures.gz and ruptures_segments.gz files for both UCERF3 branches and saves them with a .pkl extension. The ucla_plha code checks for the presence of .pkl files in the source_models directories, and uses them if they exist. Otherwise, the code uses the .gz versions of the files. 

## get_source_data(source_type, source_model, p_xyz, dist_cutoff, gmms) 

This function takes inputs of the type of source data (fault or point source), the source model (branch dependent), the site of interest, distance cutoff (e.g. 200 km), and the ground motion models of interest. For fault and point source data it returns filtered values (given distance cutoff) for magnitude, fault type (dependent on rake),rate, rjb (Joyner-Boore distance in km), rrup (rupture distance in km), rx (the horizontal distance from the top edge of the rupture in km), rx1 (in km), ry0 (in km), dip (perpendicular to strike, in radians), ztor (depth to top of rupture in km), zbor (depth to bottom of rupture in km)

Returns magnitude, fault type, rate, distance, and fault geometry terms.
    
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


## get_ground_motion_data(gmm, vs30, fault_type, rjb, rrup, rx, rx1, ry0, m, ztor, zbor, dip, z1p0, z2p5, measured_vs30)

Computes arrays containing mean and standard deviation of the natural log of a ground
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

## get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config)

Computes log-normal cumulative distribution functions for factor of safety of liquefaction

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
        config (dict): A Python dictionary read from the case insensitive config file
    
    Returns:
        fsl_cdfs (Numpy ndarray, dtype=float): Numpy array of cumulative distribution functions representing
            down-crossing rate of fsl for each event, length = N x L
        eps (Numpy ndarray, dtype=float): Numpy array of epsilon values, representing number of standard deviations
            of the natural log of fsl relative to the mean of the natural log of fsl, length = N x L
        
    Notes:
        N = number of events
        L = number of fsl values at which to compute hazard

## get_disagg(hazards, m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges)

This function calculates the seismic and liquefaction hazard disaggregation  (analyze how different earthquake scenarios (magnitude, distance, etc.) contribute to the overall seismic and liquefaction hazard at a specific location). It takes Numpy arrays of hazard values(size of number of IM values x number of events), magnitude values, distance "r" values, epsilon values, and lastly magnitude, distance, and epsilon bin edges. It returns disaggregation values for the desired bins.

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

## get_hazard(config_file)
Reads config file and runs PSHA and PLHA

    Inputs:
        config_file (string): Filename, including path, of config file. Must follow the schema
            defined by ucla_plha_schema.json. See documentation for more thorough documentation
            of the config file. The config file is case insensitive. 
    
    Returns:
        output (dict): Python dictionary containing output of analysis. The output contains all
            of the inputs for preservation, along with the hazard curve(s) and any requested
            disaggregation data. See documentation for more thorough description of output.

## References

Abrahamson, N. A., Silva, W. J., and Kamai, R. (2014). Summary of the ASK14 ground‐motion relation for active crustal regions, <em>Earthquake Spectra</em>, 30(3), 1025-1055, DOI: [10.1193/070913EQS198M](https://doi.org/10.1193/070913EQS198M)

Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014). NGA‐West2 equations for predicting PGA, PGV, and 5%-damped PSA for shallow crustal earthquakes,
<em>Earthquake Spectra</em>, 30(3), 1057--1085, DOI: [10.1193/070113EQS184M](https://doi.org/10.1193/070113EQS184M)

Boulanger, R. W. and Idriss, I. M. (2012). Probabilistic Standard Penetration Test–Based Liquefaction–Triggering Procedure. <em>J. Geotech. Geoenviron. Eng.</em>, 138(10):1185–1195. DOI: [10.1061/(ASCE)GT.1943-5606.0000700](https://doi.org/10.1061/(ASCE)GT.1943-5606.0000700)

Boulanger, R. W. and Idriss, I. M. (2016). CPT-Based Liquefaction Triggering Procedure, <em>Journal of Geotechnical and Geoenvironmental Engineering</em>, 142(2), 04015065,
DOI: [10.1061/(ASCE)GT.1943-5606.0001388](https://doi.org/10.1061/(ASCE)GT.1943-5606.0001388)

Campbell, K. W. and Bozorgnia, Y. (2014). NGA‐West2 ground motion model for the average horizontal components of PGA, PGV, and 5%-damped linear acceleration response spectra,
<em>Earthquake Spectra</em>, 30(3), 1087-1115, DOI: [10.1193/062913EQS175M](https://doi.org/10.1193/062913EQS175M)

Cetin, K.O., Seed, R.B., Kayen, R.E., Moss, R.E.S., Bilge, H.T., Ilgac, M., Chowdhury, K., et al. (2018). SPT‐based probabilistic and deterministic assessment of seismic soil liquefaction triggering hazard. <em>Soil Dynamics and Earthquake Engineering</em>, 115(12), 698-709. DOI: [10.1016/j.soildyn.2018.09.012](https://doi.org/10.1016/j.soildyn.2018.09.012)

Chiou, B. S.‐J. and Youngs, R. R. (2014). Update of the Chiou and Youngs NGA model for the average horizontal component of peak ground motion and response spectra,
<em>Earthquake Spectra</em>, 30(3), 1117-1153, DOI: [10.1193/072813EQS219M](https://doi.org/10.1193/072813EQS219M)

Field, E.H., Biasi, G.P., Bird, P., Dawson, T.E., Felzer, K.R., Jackson, D.D., Johnson, K.M., Jordan, T.H., Madden, C., Michael, A.J., Milner, K.R., Page, M.T., Parsons, T., Powers, P.M., Shaw, B.E., Thatcher, W.R., Weldon, R.J.II, and Zeng, Y. (2013). Uniform California earthquake rupture forecast, version 3 (UCERF3)—The time-independent model,
<em>U.S. Geological Survey Open-File Report</em>, 2013-1165, DOI: [10.3133/ofr20131165](https://doi.org/10.3133/ofr20131165)

Moss, R.E.S., Seed, R.B., Kayen, R.E., Stewart, J.P., Der Kiureghian, A., and Cetin, K.O. (2006). CPT-Based Probabilistic and Deterministic Assessment of In Situ Seismic Soil Liquefaction Potential, <em>Journal of Geotechnical and Geoenvironmental Engineering</em> 132(8), 1032-1051. DOI: [10.1061/(ASCE)1090-0241(2006)132:8(1032)](https://doi.org/10.1061/(ASCE)1090-0241(2006)132:8(1032))

Ulmer, K.J., Hudson, K.S., Brandenberg, S.J., Zimmaro, P., Pretell, R., Carlton, J.B., Kramer, S.L., and Stewart, J.P. (2024). Next Generation Liquefaction models for susceptibility, triggering, and manifestation, Rev. 1, <em>Research Information Letter, Office of Nuclear Regulatory Research, United States Nuclear Regulatory Commission</em>, Report No. RIL2024-13, [link](https://www.nrc.gov/docs/ML2435/ML24353A158.pdf)
