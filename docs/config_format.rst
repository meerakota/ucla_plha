Config File Format
==================

**ucla_plha** uses a `Javascript Object Notation (JSON)`_ input file to configure
the analysis, and uses `JSON Schema`_ to validate user inputs. Example config files
are provided in the :doc:`examples` section of this documentation, and the schema is
defined by the `ucla_plha_schema.json`_ file in the ``src/ucla_plha`` directory of
the GitHub repository. We validate that the data types are correct and that all
required fields have been specified. 

We thought about using a more human-readable input file format such as YAML, but
decided to use JSON instead since our code is intended for HPC applications. Both
JSON and YAML are human readable, but JSON is better for machine-parsing, and we
prioritized that in our code.

.. _JavaScript Object Notation (JSON): https://www.json.org/json-en.html
.. _JSON Schema: https://json-schema.org/
.. _PEP 8: https://peps.python.org/pep-0008/
.. _ucla_plha_schema.json: https://github.com/meerakota/ucla_plha/blob/main/src/ucla_plha/ucla_plha_schema.json


Top-level keys
--------------

The top-level keys are all JSON objects. If the ``"liquefaction_models"`` key and
its associated JSON objects are excluded from the config file, the code will perform
only a PSHA without performing a PLHA.

.. list-table::
   :header-rows: 1
   :widths: 30 55 15

   * - key
     - description
     - required
   * - ``"site"``
     - site location, shear-wave velocity, and optional depth parameters
     - yes
   * - ``"constraints"``
     - distance and magnitude limits for the analysis
     - no
   * - ``"source_models"``
     - earthquake source models (fault and point sources)
     - yes
   * - ``"ground_motion_models"``
     - ground motion models and their weights
     - yes
   * - ``"liquefaction_models"``
     - liquefaction triggering models and soil inputs; omit to run PSHA only
     - no
   * - ``"profiles"``
     - several soil profiles, each with its liquefaction models (instead of
       ``"liquefaction_models"``)
     - no
   * - ``"output"``
     - requested hazard outputs and optional disaggregation
     - yes


``"site"`` keys
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 45 18 17

   * - key
     - definition
     - units
     - required
   * - ``"latitude"``
     - latitude
     - decimal degrees
     - yes
   * - ``"longitude"``
     - longitude
     - decimal degrees
     - yes
   * - ``"elevation"``
     - elevation
     - km
     - yes
   * - ``"vs30"``
     - time-averaged shear-wave velocity in the upper 30 m
     - m/s
     - yes
   * - ``"measured_vs30"``
     - boolean indicating whether ``vs30`` was measured (``true``) or inferred (``false``)
     - \-
     - no; default = ``false``
   * - ``"z1p0"``
     - depth to the isosurface where :math:`V_S` = 1.0 km/s; used by ASK14 and CY14
     - km
     - no; default = centered value based on :math:`V_{S30}` (i.e., :math:`\delta z_{1.0} = 0`)
   * - ``"z2p5"``
     - depth to the isosurface where :math:`V_S` = 2.5 km/s; used by CB14
     - km
     - no; default = centered value based on :math:`V_{S30}` (i.e., :math:`\delta z_{2.5} = 0`)

.. note::

   :math:`z_{1.0}` is used in the BSSA14 model but has no influence on PGA, so it
   is not accepted as an input for ``"bssa14"``.


``"constraints"`` keys
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 45 18 17

   * - key
     - definition
     - units
     - required
   * - ``"dist_cutoff"``
     - maximum source-to-site distance to include in the PSHA
     - km
     - no
   * - ``"m_min"``
     - minimum magnitude to include in the PSHA
     - \-
     - no


``"source_models"`` keys
-------------------------

The second nesting level contains two keys: ``"fault_source_models"``
and ``"point_source_models"``. The third nesting level contains the specific source
model keys listed below. The fourth nesting level contains a single key
``"weight"``.

Weights will be normalized within ``"fault_source_models"`` and within
``"point_source_models"`` such that they sum to unity.

**fault_source_models**

.. list-table::
   :header-rows: 1
   :widths: 32 51 17

   * - key
     - description
     - required
   * - ``"ucerf3_fm31"``
     - UCERF3 Fault Model 3.1
     - no
   * - ``"ucerf3_fm32"``
     - UCERF3 Fault Model 3.2
     - no

**point_source_models**

.. list-table::
   :header-rows: 1
   :widths: 40 43 17

   * - key
     - description
     - required
   * - ``"ucerf3_fm31_grid_sub_seis"``
     - UCERF3 FM3.1 gridded sub-seismogenic sources
     - no
   * - ``"ucerf3_fm31_grid_unassociated"``
     - UCERF3 FM3.1 gridded unassociated sources
     - no
   * - ``"ucerf3_fm32_grid_sub_seis"``
     - UCERF3 FM3.2 gridded sub-seismogenic sources
     - no
   * - ``"ucerf3_fm32_grid_unassociated"``
     - UCERF3 FM3.2 gridded unassociated sources
     - no


``"ground_motion_models"`` keys
--------------------------------

The second nesting level defines the ground motion model and must be one or more of
``"ask14"``, ``"bssa14"``, ``"cb14"``, ``"cy14"``. The third nesting level defines
the model parameters listed below. Weights within ``"ground_motion_models"`` will be
normalized to sum to unity.

.. list-table::
   :header-rows: 1
   :widths: 20 55 25

   * - key
     - definition
     - notes
   * - ``"weight"``
     - weight assigned to this ground motion model
     - required for all models


Distance and other options
^^^^^^^^^^^^^^^^^^^^^^^^^^

``"constraints"`` accepts ``"dist_cutoff"`` (km; applied to Rjb if a ground motion model
uses Rjb, otherwise to Rrup), ``"m_min"``, ``"truncation_level"`` (upper truncation of the
ground motion distribution in standard deviations; the USGS NSHM uses 3; default none),
``"rupture_rx"``: ``"gc2"`` (default) takes Rx and Ry0 of a fault rupture from the GC2
coordinates of Spudich and Chiou (2015) of the whole rupture, as in the NGA-West2 and
NGA-West3 flatfiles (Rx1 = Rx minus the horizontal width of the segment with the smallest
Rrup), ``"closest_section"`` takes Rx and Rx1 from the segment with the smallest Rrup, as
nshmp-lib does, and ``"minimum"`` takes the minimum over the segments (earlier versions;
because Rx is signed, this is the most footwall-side segment of the rupture); with the last
two, Ry0 is the minimum over the segments. ``"triangle_frame"``: ``"local"`` (default)
computes the Rjb and Rrup of the fault triangles in a projection on a sphere centred on the
site, so that they are measured along the Earth's surface and Rrup >= Rjb, and
``"geocentric"`` uses straight lines in Earth-centred Cartesian coordinates (earlier
versions; Rrup < Rjb for most ruptures beyond about 100 km). ``"output": {"psha":
{"source_model_hazard": true}}`` adds the hazard curve of each source model.


``"liquefaction_models"`` keys
-------------------------------

The second nesting level defines the liquefaction triggering model and must be one
or more of ``"boulanger_idriss_2012"``, ``"boulanger_idriss_2016"``,
``"cetin_et_al_2018"``, ``"moss_et_al_2006"``, ``"ngl_smt_2024"``. The third
nesting level defines the model inputs listed below.

.. list-table::
   :header-rows: 1
   :widths: 16 40 10 34

   * - key
     - definition
     - units
     - required by
   * - ``"n160"``
     - overburden- and energy-corrected SPT blow count
     - \-
     - ``"boulanger_idriss_2012"``, ``"cetin_et_al_2018"``
   * - ``"fc"``
     - fines content
     - %
     - ``"boulanger_idriss_2012"``, ``"cetin_et_al_2018"``
   * - ``"vs12"``
     - time-averaged shear-wave velocity in the upper 12 m
     - m/s
     - ``"cetin_et_al_2018"``
   * - ``"parameter_uncertainty"``
     - include the input parameter uncertainty of Cetin et al. (2018b) (default ``false``)
     - \-
     - optional, ``"cetin_et_al_2018"``
   * - ``"sigma_n160"``
     - standard deviation of ``"n160"`` (default 0)
     - \-
     - optional, ``"cetin_et_al_2018"`` with ``"parameter_uncertainty"``
   * - ``"sigma_fc"``
     - standard deviation of ``"fc"`` (default 0)
     - %
     - optional, ``"cetin_et_al_2018"`` with ``"parameter_uncertainty"``
   * - ``"cov_sigmavp"``
     - coefficient of variation of ``"sigmavp"`` (default 0)
     - \-
     - optional, ``"cetin_et_al_2018"`` with ``"parameter_uncertainty"``
   * - ``"cov_csr"``
     - coefficient of variation of the cyclic stress ratio for a given ground motion (default 0)
     - \-
     - optional, ``"cetin_et_al_2018"`` with ``"parameter_uncertainty"``
   * - ``"qc1Ncs"``
     - overburden- and fines-corrected dimensionless cone tip resistance
     - \-
     - ``"boulanger_idriss_2016"``, ``"ngl_smt_2024"``
   * - ``"qc"``
     - cone tip resistance
     - MPa
     - ``"moss_et_al_2006"``
   * - ``"fs"``
     - cone sleeve friction
     - MPa
     - ``"moss_et_al_2006"``
   * - ``"rd_amax_cap"``
     - upper limit of :math:`a_{max}` in :math:`r_d` with ``"pga_uncertainty": "conditional"`` (default 1.0; ``null``: no limit)
     - g
     - optional, ``"moss_et_al_2006"``
   * - ``"Ic"``
     - soil behavior type index
     - \-
     - ``"ngl_smt_2024"``
   * - ``"Ksat"``
     - saturation factor applied to probability of triggering
     - \-
     - ``"ngl_smt_2024"``
   * - ``"pga_uncertainty"``
     - treatment of the uncertainty of the peak acceleration (see below): ``"shared"`` (default) or ``"per_layer"`` for ``"ngl_smt_2024"``; ``"conditional"`` (default) or ``"closed_form"`` for ``"boulanger_idriss_2016"`` and ``"moss_et_al_2006"``
     - \-
     - optional, ``"ngl_smt_2024"``, ``"boulanger_idriss_2016"``, ``"moss_et_al_2006"``
   * - ``"pga_integration"``
     - integration over the peak acceleration with ``"shared"`` or ``"conditional"``, ``"grid"`` (default) or ``"quadrature"``
     - \-
     - optional, ``"ngl_smt_2024"``, ``"boulanger_idriss_2016"``, ``"moss_et_al_2006"``
   * - ``"pga_grid"``
     - parameters of the magnitude :math:`\times` PGA grid: ``"pga_min"`` (default 1e-4), ``"pga_max"`` (10), ``"d_ln_pga"`` (0.1), ``"d_m"`` (0.05), ``"d_mu"`` (0.05), ``"d_sigma"`` (0.005)
     - g for ``"pga_min"`` and ``"pga_max"``
     - optional, ``"ngl_smt_2024"``, ``"boulanger_idriss_2016"``, ``"moss_et_al_2006"``
   * - ``"ztop"``
     - depth to top of layer (array)
     - m
     - ``"ngl_smt_2024"``
   * - ``"zbot"``
     - depth to bottom of layer (array)
     - m
     - ``"ngl_smt_2024"``
   * - ``"sigmav"``
     - vertical total stress at center of layer
     - kPa
     - all models; array-valued for ``"ngl_smt_2024"``, scalar-valued for other models
   * - ``"sigmavp"``
     - vertical effective stress at center of layer
     - kPa
     - all models; array-valued for ``"ngl_smt_2024"``, scalar-valued for other models
   * - ``"depth"``
     - depth to center of layer
     - m
     - ``"boulanger_idriss_2012"``, ``"boulanger_idriss_2016"``, ``"cetin_et_al_2018"``, ``"moss_et_al_2006"``
   * - ``"pa"``
     - atmospheric pressure
     - kPa
     - all models
   * - ``"weight"``
     - weight assigned to this liquefaction model
     - \-
     - all models

``"ngl_smt_2024"`` treats the peak acceleration of an earthquake as the same in all layers of the profile (``"pga_uncertainty": "shared"``, the default). The probability of profile manifestation is computed for a given magnitude :math:`M` and PGA,

.. math::

   P(FS < FS_L \mid M, PGA) = 1 - \prod_j \left[1 - P_{MT,j}\, \Phi\left(\frac{\ln FS_L - \mu_{\ln CRR,j}(M, PGA) + \ln CSR_j(M, PGA)}{0.44}\right) P_{S,j}\, K_{sat,j}\right]^{t_j / t_c}

(with the magnitude scaling factor of that PGA), and integrated over the distribution of :math:`\ln PGA` of the earthquake, truncated as in the PSHA (``"constraints": {"truncation_level": ...}``; untruncated if it is not given). The NGL SMT model itself does not treat the ground motion uncertainty. Versions 2.x and earlier included the PGA uncertainty in each layer separately (:math:`\sigma = \sqrt{0.44^2 + \sigma_{\ln PGA}^2}`, magnitude scaling factor of the median PGA, untruncated) before the product over the layers, as if each layer had an independent PGA, which gives higher probabilities (about 5% to 15% in relative terms for single earthquakes); ``"pga_uncertainty": "per_layer"`` reproduces those results exactly.

With ``"shared"``, the integration uses a grid of the rates of the ground motions of the site over magnitude and :math:`\ln PGA` (``"pga_integration": "grid"``, the default; ``ucla_plha.pga_grid``), built once for the site and shared by the liquefaction models conditioned on PGA, so that the cost of a liquefaction model is small. With the default ``"pga_grid"`` parameters (bins of 0.1 in :math:`\ln PGA` from 1e-4 g to 10 g, magnitude nodes every 0.05, and histogram nodes of the mean and standard deviation of :math:`\ln PGA` of the earthquakes every 0.05 and 0.005), the hazard curves agree with an exact integration within about 2e-4 (relative). The disaggregation uses per-earthquake probabilities from the grid nodes. ``"pga_integration": "quadrature"`` integrates each earthquake with Gauss-Legendre quadrature (slow; for verification).

``"boulanger_idriss_2016"`` and ``"moss_et_al_2006"`` are also computed conditioned on the peak acceleration (``"pga_uncertainty": "conditional"``, the default),

.. math::

   P(FS < FS_L \mid M, PGA) = \Phi\left(\frac{\ln FS_L - \mu_{\ln CRR} + \ln CSR(M, PGA)}{\sigma_{\ln CRR}}\right)

with the model uncertainty only (:math:`\sigma_{\ln CRR}` = 0.20 and 1.632 / 7.177), integrated over the ground motion distribution of the earthquake with the truncation of the PSHA, on the same magnitude :math:`\times` PGA grid as ``"ngl_smt_2024"`` (one grid per site for these models). For Boulanger and Idriss (2016), PGA enters only through CSR, so without truncation this is the closed form of earlier versions, :math:`\Phi\left((\ln FS_L - \mu_{\ln CRR} + \mu_{\ln CSR}) / \sqrt{0.20^2 + \sigma_{\ln PGA}^2}\right)`. For Moss et al. (2006), :math:`r_d` (which depends on :math:`a_{max}`) and the :math:`\Delta q_c` adjustment of the CRR (:math:`x_1 \ln CSR + x_2`, for :math:`R_f > 0.5\%`) are evaluated at the conditioned PGA instead of the median PGA; :math:`\Delta q_c` is not allowed below 0 (their Eq. 18 bounds the fines adjustment to :math:`R_f` of 0.5 to 5% as a suppression of liquefiability; :math:`x_1 \ln CSR + x_2 < 0` only for :math:`CSR < e^{-3.84} = 0.021`, outside their data); and the CSR is not divided by the duration weighting factor :math:`DWF = 17.84 M_w^{-1.43}` (with any ``"pga_uncertainty"``), because their Eq. 20 and 21 have their own magnitude term (:math:`-0.848 \ln M_w`) and are compared with the CSR of the event (the note to Eq. 21: it "produces results for the input magnitude", whereas their Figs. 6 and 10 are for :math:`M_w` 7.5, CSR*; and the CSR of the case histories in their Table 1 is :math:`0.65\, a_{max} (\sigma_v / \sigma'_v)\, r_d`, not CSR / DWF). Earlier versions divided by the DWF as well, which double counted the magnitude effect (FS higher by the DWF: about 1.4 at :math:`M_w` 6, 1.1 at 7, 0.91 at 8); the 2.x Moss et al. (2006) results, which divided CSR by the DWF, can't be reproduced with this version (a config with ``"duration_weighting"`` is rejected). ``"pga_uncertainty": "closed_form"`` gives the untruncated closed forms of earlier versions: exactly for Boulanger and Idriss (2016), and for Moss et al. (2006) with :math:`r_d` and :math:`\Delta q_c` at the median PGA but without the DWF. The :math:`r_d` of Moss et al. (2006) (from Cetin et al. 2004) decreases with :math:`a_{max}` to the lower limit 0.2 at about 2.5 to 3.5 g and jumps to the upper limit 1.0 at about 3.5 to 4 g (its denominator changes sign), so with ``"conditional"`` :math:`r_d` is evaluated with :math:`a_{max}` at most ``"rd_amax_cap"`` (default 1 g; ``null`` for no limit; :math:`r_d` above the cap is :math:`r_d` at the cap). ``"closed_form"`` (:math:`r_d` of the median PGA) has no cap. ``plha.DEFAULT_PGA_UNCERTAINTY`` holds the default treatment of each model.

``"cetin_et_al_2018"`` includes only the model uncertainty (:math:`\sigma_\varepsilon = 2.95`) by default. With ``"parameter_uncertainty": true``, the input parameter uncertainty of Cetin et al. (2018b) is added, :math:`\sigma_{tot}^2 = \theta_7 \sigma_{input}^2 + \sigma_\varepsilon^2` with :math:`\theta_7 = 0.392` and :math:`\sigma_{input}^2 = \theta_6^2 \delta_{CSR}^2 + \sigma_{N_{1,60}}^2 (1 + \theta_1 FC)^2 + \sigma_{FC}^2 (\theta_1 N_{1,60} + \theta_4)^2 + \theta_3^2 \delta_{\sigma'_v}^2` (their Eqs. 4 and 5), where :math:`\delta_{CSR}` is ``"cov_csr"``, :math:`\sigma_{N_{1,60}}` is ``"sigma_n160"``, :math:`\sigma_{FC}` is ``"sigma_fc"``, and :math:`\delta_{\sigma'_v}` is ``"cov_sigmavp"``. Because ucla_plha integrates the ground motion variability in the hazard calculation, ``"cov_csr"`` should represent only the uncertainty in estimating CSR for a given ground motion (e.g., :math:`r_d` and the stresses), not the variability of the ground motion.

.. note::

   You cannot include keys that are not required by a given model. For example,
   ``"n160"`` cannot appear in a cone penetration test model entry.


``"output"`` keys
-----------------

The ``"output"`` object specifies whether to write the PSHA and/or PLHA output to a
JSON file, and whether to include disaggregation data. Including disaggregation data
slows down the calculations, and should be included only if needed.

.. list-table::
   :header-rows: 1
   :widths: 28 45 27

   * - key
     - definition
     - notes
   * - ``"psha"``
     - include to output a PSHA hazard curve
     - optional
   * - ``"pga"``
     - array of PGA values at which to evaluate the hazard curve
     - required if ``"psha"`` is included
   * - ``"plha"``
     - include to output a PLHA hazard curve
     - optional
   * - ``"fsl"``
     - array of factor-of-safety values at which to evaluate the hazard curve
     - required if ``"plha"`` is included
   * - ``"disaggregation"``
     - include under ``"psha"`` or ``"plha"`` to output disaggregation data
     - optional
   * - ``"magnitude_bin_edges"``
     - edges of magnitude bins for disaggregation
     - required if ``"disaggregation"`` is included
   * - ``"distance_bin_edges"``
     - edges of distance bins for disaggregation
     - required if ``"disaggregation"`` is included
   * - ``"epsilon_bin_edges"``
     - edges of epsilon bins for disaggregation. Bins include their lower edge and exclude
       their upper edge; events outside the bins are not binned (use wide outer edges, e.g.
       -1e9 and 1e9, to bin all of them)
     - required if ``"disaggregation"`` is included
   * - ``"distance_metric"``
     - distance used for the distance bins: ``"default"`` (rjb if a ground motion model uses
       rjb, otherwise rrup), ``"rrup"`` (as the USGS NSHM disaggregation), or ``"rjb"``
     - optional (``"psha"`` and ``"plha"``)
   * - ``"means"``
     - ``true`` to also output ``"disaggregation_means"``: the binned contribution (%) and
       the mean magnitude, distance, and epsilon at each PGA, in total and, with
       ``"source_model_hazard"``, per source model
     - optional (``"psha"``)
   * - ``"outputfile"``
     - filename for saving output (includes a copy of the inputs)
     - optional

.. note::

   You can specify ``"plha"`` and its ``"disaggregation"`` object without including
   ``"psha"``. In this case the PSHA will still be performed internally, but its
   results will not be returned.

   Whether or not ``"outputfile"`` is specified, the code returns a Python
   dictionary containing all results.

   The input JSON object is saved alongside the output in ``"outputfile"`` so that
   inputs and outputs are always stored together. This prevents inconsistencies if
   the config file is modified after a run.


``"profiles"`` (several soil profiles at a site)
-------------------------------------------------

A config can give several soil profiles of the site with ``"profiles"`` instead of
``"liquefaction_models"``: a list of profiles, each with an optional ``"name"`` and its
``"liquefaction_models"`` (the same keys as the top-level ``"liquefaction_models"``,
including the weights, which are normalized within each profile). The source data, ground
motions, and PSHA are computed once for the site, and only the liquefaction models are
evaluated for each profile, so a site with four profiles costs much less than four runs.
The results of each profile are exactly those of a config with its liquefaction models. A
config cannot have both ``"profiles"`` and ``"liquefaction_models"``; configs without
``"profiles"`` work as before.

.. code-block:: json

   "profiles": [
      {
         "name": "boring_1",
         "liquefaction_models": {
            "moss_et_al_2006": {"qc": 5.49, "fs": 0.02617, "sigmav": 46.55,
                                "sigmavp": 43.22, "depth": 2.85, "pa": 101.325,
                                "weight": 0.5},
            "boulanger_idriss_2016": {"qc1Ncs": 71.8, "sigmav": 40.0, "sigmavp": 28.0,
                                      "depth": 2.2, "pa": 101.325, "weight": 0.5}
         }
      },
      {
         "name": "boring_2",
         "liquefaction_models": {
            "boulanger_idriss_2016": {"qc1Ncs": 95.0, "sigmav": 40.0, "sigmavp": 28.0,
                                      "depth": 2.2, "pa": 101.325, "weight": 1.0}
         }
      }
   ]

The PLHA output is then ``"output": {"profiles": [{"name": "boring_1", "plha": {...}},
{"name": "boring_2", "plha": {...}}]}``, in the order of the config, with each ``"plha"`` as
the ``"plha"`` output of a single profile (including the disaggregation if requested),
instead of ``"output": {"plha": {...}}``. Names are lowercase (the config is
case-insensitive); the default name of a profile is its index (``"0"``, ``"1"``, ...).


Running many sites
------------------

``plha.get_hazard`` also accepts a config dict (as read from a config file) instead of a file
name (``"outputfile": "default"`` then needs a file name instead). To run many sites, read and
prepare the source models once with ``plha.load_source_models`` and pass them to
``get_hazard``:

.. code-block:: python

   from ucla_plha import plha

   models = plha.load_source_models(config)  # reads and prepares the source models once
   for site in sites:
       config["site"] = site
       output = plha.get_hazard(config, models=models)

The results are identical to those of ``get_hazard`` without ``models``. The source model
arrays are read-only, so worker processes created with fork after ``load_source_models``
share them.

``python -m ucla_plha.batch`` runs a table of sites (CSV or Parquet with ``site_id``,
``latitude``, ``longitude``, ``vs30``, optional site parameters (``elevation``, ``measured_vs30``, ``z1p0``, ``z2p5``), and an optional
``profiles_file`` with the ``"profiles"`` of each site) with a base config, loads the source
models once, and processes the sites with a pool of worker processes (fork on Linux, so the
workers share the source model arrays; spawn elsewhere). Each chunk of sites is written by its
worker to one compact file (npz, or HDF5 with ``--format hdf5``) with the site ids, the PSHA
curves, and the PLHA curves of every profile; existing chunks are skipped, so an interrupted
run can be restarted:

.. code-block:: bash

   python -m ucla_plha.batch --sites sites.csv --config base.json --out-dir results \
       --processes 128 --chunk-size 100 [--start 0 --max-sites 100000]

On a cluster node (e.g. TACC Stampede3 or Vista) run one process per core
(``--processes`` = number of cores, with ``OMP_NUM_THREADS=1``), and split a large table
over a SLURM job array with ``--start`` and ``--max-sites`` (e.g.
``--start $((SLURM_ARRAY_TASK_ID * 100000)) --max-sites 100000``), so that each array task
writes its own chunk files. ``ucla_plha.batch.read_results`` reads the chunk files. mpi4py
could distribute chunks over the nodes of one job instead, but is not needed. See the
docstring of ``ucla_plha/batch.py``.

GC2 cache
^^^^^^^^^

The site-independent GC2 set-up of the fault source models (Rx and Ry0 with
``"rupture_rx": "gc2"``, the default) takes about 6 s for UCERF3 FM3.1 and 7 s for FM3.2 in
each new process, and is too large to ship with the package (147 and 177 MB). It is written to a cache directory the
first time it is computed and read by later processes (0.1 to 0.3 s per model):
``UCLA_PLHA_CACHE_DIR`` (default ``~/.cache/ucla_plha``, or ``%LOCALAPPDATA%\ucla_plha\Cache``
on Windows; ``none`` disables the cache). Run the precompute step once after installing, so that
the jobs of a cluster run do not each compute it:

.. code-block:: bash

   export UCLA_PLHA_CACHE_DIR=$WORK/ucla_plha_cache   # TACC: not in $HOME (small quota)
   python -m ucla_plha.precompute                     # both UCERF3 fault models, about 325 MB

and use the same ``UCLA_PLHA_CACHE_DIR`` in the job scripts. The files are keyed on a hash of
the model geometry files, the cache format, and the numpy version, written atomically (a
temporary file in the cache directory that is renamed), and checked when they are read (key,
array sizes, CRC-32 of the npz members); a missing, stale, or corrupt file is computed again,
and the results are identical with and without the cache (``ucla_plha.gc2_cache``). With the
cache, the first ``get_hazard`` of a new process for ``examples/config10.json`` without its
disaggregation takes about 4 s instead of 19 s.


Full config structure
---------------------

The following shows the complete structure of a config file. Note that keys are
lowercase following `PEP 8`_ convention, but the config file is case-insensitive
so ``N160`` works equally well as ``n160``. Optional inputs are shown with their
default values where applicable.

.. code-block:: json

   {
      "site": {
         "latitude": 34.0,
         "longitude": -118.0,
         "elevation": 0.0,
         "vs30": 350.0,
         "measured_vs30": false,
         "z1p0": 0.05,
         "z2p5": 0.5
      },
      "constraints": {
         "dist_cutoff": 200,
         "m_min": 5.0
      },
      "source_models": {
         "fault_source_models": {
            "ucerf3_fm31": {"weight": 0.5},
            "ucerf3_fm32": {"weight": 0.5}
         },
         "point_source_models": {
            "ucerf3_fm31_grid_sub_seis":     {"weight": 0.25},
            "ucerf3_fm31_grid_unassociated": {"weight": 0.25},
            "ucerf3_fm32_grid_sub_seis":     {"weight": 0.25},
            "ucerf3_fm32_grid_unassociated": {"weight": 0.25}
         }
      },
      "ground_motion_models": {
         "ask14":  {"weight": 0.25},
         "bssa14": {"weight": 0.25},
         "cb14":   {"weight": 0.25},
         "cy14":   {"weight": 0.25}
      },
      "liquefaction_models": {
         "boulanger_idriss_2012": {
            "n160": 10.0, "fc": 5.0,
            "sigmav": 100.0, "sigmavp": 60.0, "depth": 5.0, "pa": 101.325,
            "weight": 0.25
         },
         "boulanger_idriss_2016": {
            "qc1Ncs": 80.0,
            "sigmav": 100.0, "sigmavp": 60.0, "depth": 5.0, "pa": 101.325,
            "weight": 0.25
         },
         "cetin_et_al_2018": {
            "n160": 10.0, "fc": 5.0, "vs12": 150.0,
            "sigmav": 100.0, "sigmavp": 60.0, "depth": 5.0, "pa": 101.325,
            "weight": 0.25
         },
         "moss_et_al_2006": {
            "qc": 3.0, "fs": 0.05,
            "sigmav": 100.0, "sigmavp": 60.0, "depth": 5.0, "pa": 101.325,
            "weight": 0.25
         },
         "ngl_smt_2024": {
            "ztop":    [4.5, 5.0],
            "zbot":    [5.0, 5.5],
            "qc1Ncs":  [80.0, 85.0],
            "Ic":      [1.8, 1.9],
            "sigmav":  [100.0, 110.0],
            "sigmavp": [60.0, 65.0],
            "Ksat":    [1.0, 1.0],
            "pa": 101.325,
            "weight": 0.0
         }
      },
      "output": {
         "psha": {
            "pga": [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0],
            "disaggregation": {
               "magnitude_bin_edges": [5.0, 6.0, 7.0, 8.0],
               "distance_bin_edges":  [0, 25, 50, 100, 200],
               "epsilon_bin_edges":   [-3, -2, -1, 0, 1, 2, 3]
            }
         },
         "plha": {
            "fsl": [0.5, 0.75, 1.0, 1.25, 1.5, 2.0],
            "disaggregation": {
               "magnitude_bin_edges": [5.0, 6.0, 7.0, 8.0],
               "distance_bin_edges":  [0, 25, 50, 100, 200],
               "epsilon_bin_edges":   [-3, -2, -1, 0, 1, 2, 3]
            }
         },
         "outputfile": "output.json"
      }
   }
