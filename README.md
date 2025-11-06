# ucla_plha: A Probabilistic Liquefaction Hazard Assessment Python Package

The ucla_plha Python package performs probabilistic liquefaction hazard analysis (PLHA) to quantify the annual down-crossing rate (or alternatively, the annual rate of non-exceedance) of factor of safety of liquefaction. The code integrates probabilistic seismic hazard analysis (PSHA) within the same computational engine as PLHA. The ucla_plha_code performs PLHA inside the hazard integral, which distinguishes it from other codes that decouple the PSHA from the PLHA. Currently, the code uses the Uniform California Earthquake Rupture Forecast, Version 3 (UCERF3, Field et al. 2013) source model, and is therefore applicable for shallow crustal earthquakes in California. The code uses the four NGAWest2 ground motion models that model site response using $V_{S30}$, including Abrahamson et al. (2014), Boore et al. (2014), Campbell and Bozorgnia (2014), and Chiou and Youngs (2014). The code currently uses five probabilistic liquefaction models, including two standard penetration test models (Boulanger and Idriss 2012, and Cetin et al. 2018) and three cone penetration test models (Moss et al. 2006, Boulanger and Idriss 2016, and Ulmer et al. 2024). More thorough documentation is provided in the [wiki](https://github.com/meerakota/ucla_plha/wiki).

## Installation
> pip install ucla_plha

## Documentation
> https://github.com/meerakota/ucla_plha/wiki

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
