Limitations
===========

**ucla_plha** was developed specifically to evaluate the influence of sea level rise on liquefaction in coastal areas of California,
and therefore has a number of limitations that currently prevent its application in other regions. The code is open source, and 
can therefore in principle be modified in the future to address the limitations.


#. The code currently uses the Uniform California Earthquake Rupture Forecast, Version 3 (`UCERF3, Field et al. 2013`_) source model, 
   and is therefore applicable for shallow crustal earthquakes in California.  
#. The code currently uses the four NGAWest2 ground motion models that use VS30 in site response modeling(`Abrahamson et al. 2014`_, `Boore et al. 2014`_, 
   `Campbell and Bozorgnia 2014`_, and `Chiou and Youngs 2014`_), and is therefore applicable only to shallow crustal earthquakes.  
#. The code currently supports peak horizontal acceleration, and is therefore usable only for liquefaction models that use PGA.  
#. The code currently uses five probabilistic liquefaction models, including two standard penetration test models 
   (`Boulanger and Idriss 2012`_, and `Cetin et al. 2018`_) and three cone penetration test models 
   (`Moss et al. 2006`_, `Boulanger and Idriss 2016`_, and `Ulmer et al. 2024`_).  


.. _UCERF3, Field et al. 2013: https://wgcep.org/UCERF3.html
.. _Abrahamson et al. 2014: https://doi.org/10.1193/070913EQS198M
.. _Boore et al. 2014: https://doi.org/10.1193/070113EQS184M
.. _Campbell and Bozorgnia 2014: https://doi.org/10.1193/062913EQS175M
.. _Chiou and Youngs 2014: https://doi.org/10.1193/072813EQS219M
.. _Boulanger and Idriss 2012: https://doi.org/10.1061/(ASCE)GT.1943-5606.0000700
.. _Cetin et al. 2018: https://doi.org/10.1016/j.soildyn.2018.09.012
.. _Moss et al. 2006: https://doi.org/10.1061/(ASCE)1090-0241(2006)132:8(1032)
.. _Boulanger and Idriss 2016: https://doi.org/10.1061/(ASCE)GT.1943-5606.0001388
.. _Ulmer et al. 2024: https://www.nrc.gov/docs/ML2435/ML24353A158.pdf
