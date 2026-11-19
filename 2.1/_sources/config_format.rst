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
   * - ``"Ic"``
     - soil behavior type index
     - \-
     - ``"ngl_smt_2024"``
   * - ``"Ksat"``
     - saturation factor applied to probability of triggering
     - \-
     - ``"ngl_smt_2024"``
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
     - edges of epsilon bins for disaggregation
     - required if ``"disaggregation"`` is included
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
