Functions
=========

The functions available in the **ucla_plha** package are summarized below, followed
by a detailed description of each. Some general notes on the functions:
- :math:`N` = number of events
- :math:`L` = number of FSL or intensity measure values at which to compute hazard curves
- :math:`M` = number of magnitude bins for disaggregation
- :math:`R` = number of distance bins for disaggregation
- :math:`E` = number of epsilon bins for disaggregation

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - function
     - description
   * - :func:`decompress_ucerf3_source_data`
     - Decompress the UCERF3 source-model files
   * - :func:`get_source_data`
     - Read source-model files and compute rupture distances and geometry
   * - :func:`get_ground_motion_data`
     - Evaluate the four :math:`V_{S30}`-dependent ground motion models
   * - :func:`get_liquefaction_cdfs`
     - Evaluate probabilistic liquefaction triggering models
   * - :func:`get_disagg`
     - Compute seismic and liquefaction hazard disaggregation
   * - :func:`get_hazard`
     - Run the full PSHA and PLHA from a config file

----

.. py:function:: decompress_ucerf3_source_data()

   Some of the input UCERF3 source-model files are quite large and are compressed as
   NumPy archives to keep the package installation efficient. Within the
   ``ucerf3_fm31`` and ``ucerf3_fm32`` directories under
   ``source_models/fault_source_models``, the ``ruptures.npz`` and
   ``ruptures_segments.npz`` files are compressed as (``.npz``) NumPy archives. An .npz file
   is a NumPy specific zip archive format used to save multiple NumPy arrays into a single file on disk. Loading
   the compressed archives at run time is slower than loading uncompressed data, so
   this function reads the ``.npz`` archives for both UCERF3 branches and writes
   uncompressed ``.npy`` versions. The **ucla_plha** code checks for the ``.npy``
   files first and uses them if present, otherwise falling back to the ``.npz``
   archives. It only needs to be run once.


.. py:function:: get_source_data(source_type, source_model, p_xyz, dist_cutoff, gmms)

   This function takes inputs of the type of source data (fault or point source),
   the source model (branch dependent), the site of interest, distance cutoff (e.g. 200 km),
   and the ground motion models of interest. For fault and point source data it returns
   filtered values (given distance cutoff) for magnitude, fault type (dependent on rake),
   rate, :math:`R_{JB}` (Joyner-Boore distance in km), :math:`R_{rup}` (rupture distance
   in km), :math:`R_x` (the horizontal distance from the top edge of the rupture in km),
   :math:`R_{x1}` (in km), :math:`R_{y0}` (in km), dip (perpendicular to strike, in degrees),
   :math:`Z_{tor}` (depth to top of rupture in km), and :math:`Z_{bor}` (depth to bottom of
   rupture in km).

   Returns magnitude, fault type, rate, distance, and fault geometry terms.

   :param str source_type: Either ``"fault_source_models"`` or ``"point_source_models"``.
   :param str source_model: Directory for the source model within the source-type
      directory. Currently either ``"ucerf3_fm31"`` or ``"ucerf3_fm32"``.
   :param p_xyz: Array containing the x, y, z coordinates of the point of interest,
      length 3.
   :type p_xyz: numpy.ndarray of float
   :param float dist_cutoff: Maximum distance to consider in the seismic hazard analysis.
   :param gmms: Ground motion models to use, one or more of ``"ask14"``, ``"bssa14"``,
      ``"cb14"``, ``"cy14"``.
   :type gmms: array of str

   :returns: A tuple of arrays, each of length :math:`N`:

      * **m** -- magnitudes.
      * **fault_type** -- fault types (1 = reverse, 2 = normal, 3 = strike slip).
      * **rate** -- annual rate of occurrence of each event.
      * **rjb** -- :math:`R_{JB}`, Joyner-Boore distance between the site and each rupture, in km.
      * **rrup** -- :math:`R_{rup}`, rupture distance between the site and each rupture, in km.
      * **rx** -- :math:`R_x`, distance to the surface projection of the top of each rupture,
        measured perpendicular to strike, in km.
      * **rx1** -- :math:`R_{x1}`, distance to the surface projection of the bottom of each rupture,
        measured perpendicular to strike, in km.
      * **ry0** -- :math:`R_{y0}`, distance to the surface projection of each rupture, measured
        parallel to strike, in km.
      * **dip** -- dip angle of each rupture, in degrees.
      * **ztor** -- :math:`Z_{tor}`, depth to the top of each rupture, in km.
      * **zbor** -- :math:`Z_{bor}`, depth to the bottom of each rupture, in km.


.. py:function:: get_ground_motion_data(gmm, vs30, fault_type, rjb, rrup, rx, rx1, ry0, m, ztor, zbor, dip, z1p0, z2p5, measured_vs30)

   Computes arrays containing the mean and standard deviation of the natural log of a
   ground motion intensity measure.

   :param str gmm: Ground motion model. One of ``"ask14"``, ``"bssa14"``, ``"cb14"``, ``"cy14"``.
   :param float vs30: :math:`V_{S30}`, time-averaged shear-wave velocity in the upper 30 m, in m/s.
   :param bool measured_vs30: Whether :math:`V_{S30}` is measured (``True``) or inferred (``False``).
   :param float z1p0: :math:`Z_{1.0}`, depth to the 1.0 km/s shear-wave isosurface, in km.
   :param float z2p5: :math:`Z_{2.5}`, depth to the 2.5 km/s shear-wave isosurface, in km.
   :param fault_type: Fault type for each rupture (1 = reverse, 2 = normal, 3 = strike slip), length :math:`N`.
   :param rjb: :math:`R_{JB}`, Joyner-Boore distances, length :math:`N`.
   :param rrup: :math:`R_{rup}`, rupture distances, length :math:`N`.
   :param rx: :math:`R_x`, distance to the surface projection of the top of each rupture, perpendicular to strike, length :math:`N`.
   :param rx1: :math:`R_{x1}`, distance to the surface projection of the bottom of each rupture, perpendicular to strike, length :math:`N`.
   :param ry0: :math:`R_{y0}`, distance to the surface projection of each rupture, parallel to strike, length :math:`N`.
   :param m: Magnitudes, length :math:`N`.
   :param ztor: :math:`Z_{tor}`, depth to the top of each rupture, in km, length :math:`N`.
   :param zbor: :math:`Z_{bor}`, depth to the bottom of each rupture, in km, length :math:`N`.
   :param dip: Dip angle of each rupture, in degrees, length :math:`N`.

   :returns:
      * **mu_ln_pga** -- :math:`\mu_{\ln PGA}`, mean of the natural log of the ground motion intensity measure, length :math:`N`.
      * **sigma_ln_pga** -- :math:`\sigma_{\ln PGA}`, standard deviation of the natural log of the intensity measure, length :math:`N`.


.. py:function:: get_liquefaction_cdfs(m, mu_ln_pga, sigma_ln_pga, fsl, liquefaction_model, config)

   Computes log-normal cumulative distribution functions for the factor of safety of
   liquefaction.

   :param m: Magnitudes, length :math:`N`.
   :param mu_ln_pga: :math:`\mu_{\ln PGA}`, mean of the natural log of the ground motion intensity measure, length :math:`N`.
   :param sigma_ln_pga: :math:`\sigma_{\ln PGA}`, standard deviation of the natural log of the intensity measure, length :math:`N`.
   :param fsl: Factor-of-safety values at which to compute the liquefaction hazard curve, length :math:`L`.
   :param str liquefaction_model: One of ``"cetin_et_al_2018"``, ``"moss_et_al_2006"``,
      ``"boulanger_idriss_2016"``, ``"boulanger_idriss_2012"``, ``"ngl_smt_2024"``.
   :param dict config: Dictionary read from the (case-insensitive) config file.

   :returns:
      * **fsl_cdfs** -- cumulative distribution functions representing the down-crossing
        rate of :math:`FS_L` for each event, shape :math:`N \times L`.
      * **eps** -- :math:`\varepsilon` values (number of standard deviations of :math:`\ln(FS_L)`
        relative to its mean), shape :math:`N \times L`.


.. py:function:: get_disagg(hazards, m, r, eps, m_bin_edges, r_bin_edges, eps_bin_edges)

   This function calculates the seismic and liquefaction hazard disaggregation (i.e., how
   different earthquake scenarios (magnitude, distance, etc.) contribute to the overall
   seismic and liquefaction hazard at a specific location).
   It takes NumPy arrays of hazard values (size of number of IM values x number of events),
   magnitude values, distance values, :math:`\varepsilon` values, and lastly magnitude, distance,
   and :math:`\varepsilon` bin edges. It returns disaggregation values for the desired bins.

   :param hazards: Hazard values, shape :math:`L \times N`.
   :param m: Magnitude values, length :math:`N`.
   :param r: Distance values, length :math:`N`.
   :param eps: :math:`\varepsilon` values, shape :math:`L \times N`.
   :param m_bin_edges: Magnitude bin edges, length :math:`M + 1`.
   :param r_bin_edges: Distance bin edges, length :math:`R + 1`.
   :param eps_bin_edges: :math:`\varepsilon` bin edges, length :math:`E + 1`.

   :returns: **disagg** -- contribution to the hazard within each magnitude, distance,
      and :math:`\varepsilon` bin, for each PGA (PSHA) or :math:`FS_L` (PLHA) value, shape
      :math:`L \times M \times R \times E`.


.. py:function:: get_hazard(config_file)

   Read a config file and run the PSHA and/or PLHA. This is the main entry point of the
   package.

   :param str config_file: Filename, including path, of the config file. Must follow
      the schema defined by ``ucla_plha_schema.json`` (see :doc:`config_format`). The
      config file is case-insensitive.

   :returns: **output** -- dictionary containing the analysis output. It preserves all
      of the inputs and adds the hazard curve(s) and any requested disaggregation data.


References
----------

.. [ASK14] Abrahamson, N. A., Silva, W. J., and Kamai, R. (2014). Summary of the ASK14 ground-motion relation for active crustal regions, *Earthquake Spectra*, 30(3), 1025-1055. DOI: `10.1193/070913EQS198M <https://doi.org/10.1193/070913EQS198M>`_

.. [BSSA14] Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014). NGA-West2 equations for predicting PGA, PGV, and 5%-damped PSA for shallow crustal earthquakes, *Earthquake Spectra*, 30(3), 1057-1085. DOI: `10.1193/070113EQS184M <https://doi.org/10.1193/070113EQS184M>`_

.. [BI12] Boulanger, R. W. and Idriss, I. M. (2012). Probabilistic Standard Penetration Test-Based Liquefaction-Triggering Procedure, *Journal of Geotechnical and Geoenvironmental Engineering*, 138(10), 1185-1195. DOI: `10.1061/(ASCE)GT.1943-5606.0000700 <https://doi.org/10.1061/(ASCE)GT.1943-5606.0000700>`_

.. [BI16] Boulanger, R. W. and Idriss, I. M. (2016). CPT-Based Liquefaction Triggering Procedure, *Journal of Geotechnical and Geoenvironmental Engineering*, 142(2), 04015065. DOI: `10.1061/(ASCE)GT.1943-5606.0001388 <https://doi.org/10.1061/(ASCE)GT.1943-5606.0001388>`_

.. [CB14] Campbell, K. W. and Bozorgnia, Y. (2014). NGA-West2 ground motion model for the average horizontal components of PGA, PGV, and 5%-damped linear acceleration response spectra, *Earthquake Spectra*, 30(3), 1087-1115. DOI: `10.1193/062913EQS175M <https://doi.org/10.1193/062913EQS175M>`_

.. [CEA18] Cetin, K. O., Seed, R. B., Kayen, R. E., Moss, R. E. S., Bilge, H. T., Ilgac, M., and Chowdhury, K. (2018). SPT-based probabilistic and deterministic assessment of seismic soil liquefaction triggering hazard, *Soil Dynamics and Earthquake Engineering*, 115, 698-709. DOI: `10.1016/j.soildyn.2018.09.012 <https://doi.org/10.1016/j.soildyn.2018.09.012>`_

.. [CY14] Chiou, B. S.-J. and Youngs, R. R. (2014). Update of the Chiou and Youngs NGA model for the average horizontal component of peak ground motion and response spectra, *Earthquake Spectra*, 30(3), 1117-1153. DOI: `10.1193/072813EQS219M <https://doi.org/10.1193/072813EQS219M>`_

.. [UCERF3] Field, E. H., Biasi, G. P., Bird, P., Dawson, T. E., Felzer, K. R., Jackson, D. D., Johnson, K. M., Jordan, T. H., Madden, C., Michael, A. J., Milner, K. R., Page, M. T., Parsons, T., Powers, P. M., Shaw, B. E., Thatcher, W. R., Weldon, R. J. II, and Zeng, Y. (2013). Uniform California earthquake rupture forecast, version 3 (UCERF3) — The time-independent model, *U.S. Geological Survey Open-File Report*, 2013-1165. DOI: `10.3133/ofr20131165 <https://doi.org/10.3133/ofr20131165>`_

.. [MEA06] Moss, R. E. S., Seed, R. B., Kayen, R. E., Stewart, J. P., Der Kiureghian, A., and Cetin, K. O. (2006). CPT-Based probabilistic and deterministic assessment of in situ seismic soil liquefaction potential, *Journal of Geotechnical and Geoenvironmental Engineering*, 132(8), 1032-1051. DOI: `10.1061/(ASCE)1090-0241(2006)132:8(1032) <https://doi.org/10.1061/(ASCE)1090-0241(2006)132:8(1032)>`_

.. [NGL24] Ulmer, K. J., Hudson, K. S., Brandenberg, S. J., Zimmaro, P., Pretell, R., Carlton, J. B., Kramer, S. L., and Stewart, J. P. (2024). Next Generation Liquefaction models for susceptibility, triggering, and manifestation, Rev. 1, *Research Information Letter, Office of Nuclear Regulatory Research, United States Nuclear Regulatory Commission*, Report No. RIL2024-13. `link <https://www.nrc.gov/docs/ML2435/ML24353A158.pdf>`_
