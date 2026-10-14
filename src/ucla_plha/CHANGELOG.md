Unreleased
* cetin_et_al_2018: optional input parameter uncertainty of Cetin et al. (2018b, MethodsX 5, 1556-1575; Eqs. 4 and 5). It is off by default ("parameter_uncertainty": false, the earlier behavior); with "parameter_uncertainty": true, "sigma_n160", "sigma_fc", "cov_sigmavp", and "cov_csr" (default 0) are combined with the model uncertainty as sigma_tot^2 = 0.392 sigma_input^2 + 2.95^2.

Version 2.0.0 (November 6, 2025)
* Completed documentation
* Updated ngl_smt_2024 model to accept arrays for ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, and Ksat instead of reading them from a CSV file.
* Changed idriss_boulanger_2012 to boulanger_idriss_2012
* Updated ucla_plha_schema.json file to reflect changes to schema
* Updated example config.json file to reflect changes to schema