Coupled vs. Decoupled approaches
================================
**ucla_plha** uses a coupled approach, which means that the PLHA calculations are performed inside the hazard integral, as indicated
by Eq. 1. An alternative method is to decouple the PSHA from the PLHA by first performing a PSHA and disaggregation, and subsequently
numerically convolving the seismic hazard curve with the liquefaction model. This is equivalent to pulling the liquefaction integral
out of the hazard integration, as shown by Eq. 2., which is generally simplified as represented by Eq. 3

A list of pros and cons of each approach is provided below so users can decide which approach is best suited to their application.

Coupled pros / cons
-------------------

The pros of the coupled approach are 

#.  it is computationally efficient when ground motion demand and liquefaction capacity are represented by probability density 
    functions of the same functional form (i.e., log-normal) because convolution can be performed quickly in closed-form 
#.  it inherently accounts for the influence of magnitude on liquefaction because liquefaction assessment is performed for every 
    event considered in the PSHA 
#.  it is mathematically rigorous because liquefaction fragility depends on magnitude, and cannot be pulled out of the integral as
    is done in Eq. 2.

The cons of the coupled approach are 

#.  computer code must be developed to integrate PSHA with PLHA (which is the purpose of the **ucla_plha** 
    package)
#.  convolution becomes computationally expensive when ground motion demand and liquefaction capacity for a given 
    event have different distribution functions. 

Decoupled pros / cons
---------------------

The pros of the decoupled approach are 

#.  the PSHA can be performed using any available seismic hazard code, as long as a disaggregation is provided
#.  the computation time is the same whether the ground motion and liquefaction fragility functions have the same functional form

The primary con of the decoupled approach is

#.  it is not mathematically rigorous because liquefaction capacity depends on magnitude, and disaggregation at multiple 
    return periods must be used to identify a weighted set of magnitudes for performing convolution. This step is time 
    consuming, involves judgment, and is a numerical approximation that asymptotically approaches the exact solution as 
    the number of considered magnitudes approaches the number of distinct events considered in the PSHA.

Equations
---------

.. math::
    :label: Eq. 1

    \Lambda\left(FS_L < fs_l\right) = \sum\limits_{i=1}^{n_{sources}} \lambda \left(\textbf{M}_i > \textbf{m}_{min}\right)
    \int\limits_{m_{min}}^{m_{max}} \int\limits_{0}^{r_{max}} \int\limits_{-\infty}^{+\infty} P\left(FS_L < fs_l | im,m,\chi\right) f_{IM}\left(im\right|m,r,\xi) \ f_R\left(r\right) \ f_M\left(m\right) \ dim \ dr \ dm  

.. math::
    :label: Eq. 2

    \Lambda\left(FS_L < fs_l\right) = \int\limits_{-\infty}^{+\infty} P\left(FS_L < fs_l | im, m, \chi \right) 
    \sum\limits_{i=1}^{n_{sources}} \lambda \left(M_i > m_{min}\right) \ \int\limits_{m_{min}}^{m_{max}} \int\limits_{0}^{r_{max}} f_{IM}\left(im\right | m,r, \xi)\ f_R\left(r\right) \ f_M\left(m\right) \ dr \ dm \ dim

.. math::
    :label: Eq. 3

    \Lambda\left(FS_L < fs_l\right) = \int\limits_{-\infty}^{+\infty} P\left(FS_L < fs_l | im, m \right) \ d\lambda\left(IM > im\right) dim

.. math::
    :label: Eq. 4

    \frac{d\lambda\left(IM > im\right)}{dim} = \sum\limits_{i=1}^{n_{sources}} \lambda \left(M_i > m_{min}\right) \int\limits_{m_{min}}^{m_{max}} \int\limits_{0}^{r_{max}} f_{IM}\left(im\right | m,r,\xi) \ f_R\left(r\right) \ f_M\left(m\right) \ dr \ dm 

Definition of Terms
-------------------

+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| variable                                       | definition                                                                                                |
+================================================+===========================================================================================================+
| :math:`FS_L`                                   | liquefaction factor of safety                                                                             |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\textbf{M}`                             | magnitude                                                                                                 |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+ 
| :math:`R`                                      | distance                                                                                                  |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\textbf{m}_{min}`                       | minimum magnitude                                                                                         |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`r_{max}`                                | maximum distance                                                                                          |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`n_{sources}`                            | number of earthquake sources                                                                              |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`IM`                                     | ground motion intensity measure                                                                           |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\chi`                                   | variables representing soil properties in liquefaction model                                              |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\xi`                                    | variables representing site conditions, style of faulting, etc.                                           |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\Lambda\left(FS_L < fs_l\right)`        | annual rate at which :math:`FS_L` drops below :math:`fs_l`                                                |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`\lambda \left(M_i > m_{min}\right)`     | annual rate of of earthquakes with :math:`\textbf{M} > \textbf{m}`                                        |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`P\left(FS_L < fs_l | im,m,\chi\right)`  | liquefaction fragility function representing probability of :math:`FS_L` dropping below :math:`fs_l`      |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`f_{IM}\left(im\right|\textbf{m},r,\xi)` | probability density function for :math:`IM` conditioned on :math:`\textbf{m}`, :math:`R`, and :math:`\xi` |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`f_R\left(r\right)`                      | probability density function for distance                                                                 |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| :math:`f_M\left(\textbf{m}\right)`             | probability density function for magnitude                                                                |
+------------------------------------------------+-----------------------------------------------------------------------------------------------------------+
