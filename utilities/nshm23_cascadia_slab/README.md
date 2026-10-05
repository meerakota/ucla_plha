# Subduction Intraslab Sources

In the 2008 NSHM, deep intraslab earthquakes in northern California, Oregon and Washington
were modeled as a single, 50 km, depth slice. In the 2014 NSHM, region dependent depth slices
were used. Early implementations kept the individual slices separate in the logic tree.
With the move to using variable depth spatial PDFs for 2023, the 2014/2018 model was refactored
to combine the depth slices.

The 2014/2018 Oregon grid is 50/50 Gaussian fixed kernel smoothing and fixed rate grid. These
branches have been combined in the OR rate files.

## Polygon Depths and IDs

| Grouping (Depth: id) | CA             | OR             | WA             |
|:---------------------|:--------------:|:--------------:|:--------------:|
| 2008                 | 50 km: 8200    |  50 km: 8210   | 50 km: 8220    |
| 2014 Shallow         | 39 km: 8201    | 42 km: 8211    | 42 km: 8221    |
| 2014 Middle          | 46 km: 8202    | 50 km: 8212    | 50 km: 8222    |
| 2014 Deep            | 60 km: 8203    | 60 km: 8213    | 60 km: 8223    |
| 2014/2018 refactored | variable: 8200 | variable: 8210 | variable: 8220 |

## Rupture Set IDs

Oregon and Washington intraslab MFDs are split into high and low branches, each with unique
*b*-values and corresponding *a*-values.

| State: MFD                       | Shallow | Middle | Deep | var  |
|:---------------------------------|:-------:|:------:|:----:|:----:|
| __CA:__ M 5.0-7.2 b=0.8          | 8204    | 8205   | 8206 | 8201 |
| __OR:__ M 6.5-7.2 b=0.4 (lo)     | 8214    | 8216   | 8218 | 8211 |
| __OR:__ M 7.2-7.5,8.0 b=0.8 (hi) | 8215    | 8217   | 8219 | 8212 |
| __WA:__ M 5.0-7.2 b=0.4  (lo)    | 8224    | 8226   | 8228 | 8221 |
| __WA:__ M 7.2-7.5,8.0 b=0.8 (hi) | 8225    | 8227   | 8229 | 8222 |

## Source Tree IDs

__CA:__ 8230  
__OR:__ 8231  
__WA:__ 8232  
