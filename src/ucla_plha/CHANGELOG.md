Unreleased
* Rx and Rx1 of fault ruptures with several fault segments are now those of the segment with the smallest Rrup, as in nshmp-lib (SystemRuptureSet), instead of the minimum over the segments, which, because Rx is signed, was the most footwall-side segment and often turned hanging-wall sites into footwall sites. The earlier convention is available with "constraints": {"rupture_rx": "minimum"}.
* nshm23_wus fault geometry: the aseismic reduction of the seismogenic width now removes the top of the dipping fault plane (the upper edge moves down-dip, as in nshmp-lib DefaultGriddedSurface) instead of moving the whole plane down, which displaced dipping faults toward their footwall by aseismicity * (lower depth - upper depth) / tan(dip) (2.5 to 2.9 km for the Puente Hills and Compton thrusts).

Version 2.0.0 (November 6, 2025)
* Completed documentation
* Updated ngl_smt_2024 model to accept arrays for ztop, zbot, qc1Ncs, Ic, sigmav, sigmavp, and Ksat instead of reading them from a CSV file.
* Changed idriss_boulanger_2012 to boulanger_idriss_2012
* Updated ucla_plha_schema.json file to reflect changes to schema
* Updated example config.json file to reflect changes to schema