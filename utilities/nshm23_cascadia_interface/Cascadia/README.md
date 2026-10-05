# Cascadia Subduction Zone (CSZ)

The Cascadia subduction zone model was updated in the Conterminous U.S. NSHM from 2008 to 2014,
remained the same in 2018, and was updated again in 2023. The Cascadia logic tree as a whole has
the ID: `3199`.

## Magnitude Frequency Distributions (MFDs)

The single (a.k.a. 'characteristic') magnitude frequency distributions for the different sections
and combinations thereof are based on the magnitude (area) scaling relations listed below,
in order of increasing implied magnitude, and are identified as follows:

- `M1`: Strasser et al. (2010)
- `M2`: Murotani et al. (2008)
- `M3`: Papazachos et al. (2004)
- `M4`: Shaw (c = 4.0)

Cascadia uses two different Gutenberg–Richter MFD's, one with `b=0.0`, and one with `b=1.0`,
both of which are used to float ruptures on the 'full' and models spanning the 3 southern sections.

For the 2014 & 2018 Conterminous U.S. NSHM update, Cascadia either ruptures in its entirety in a
so-called 'M9' event or in a series of smaller ~M8 events. Either are possible and the logic
trees for each are additive to give the overall Goldfinger et al. (2012) rate for the entire
subduction zone. Therefore there is no weighting between the 'full' and 'partial' branches of
the logic tree defined in Fig. 63 of Petersen et al. (2014). `source-groups` are used to define
additive branches and to supply scale factors.

- Original full and segmented files: `sub[0|1|2|3|4]_ch_[bot|mid|top].in`
- Original unsegmented files: `sub[0|1]_GRb[0|1]_[bot|mid|top].in`

## Section ID Groupings

The ID ranges of features used to define the geometry of current and past Casacadia models are
listed in the table below.

| Grouping              | ID Range             |
| ---------------------:|:-------------------- |
| 2008 Sections         | `[3100..3103]`       |
| 2014/2018 Sections    | `[3110..3142]`       |
| 2014/2018 Ruptures    | `[3150..3149]`       |
| 2023 Sections         | `[3200..3281]`       |
| 2023 Ruptures         | `[100000..100902]`   |

TODO: Most likely will change the rupture ID's for 2023. Replace in table.

## Rupture Reference (2023)

The table below lists the section IDs used for rupture sets on branches of the 2023
Cascadia logic tree. The Full and Unsegmented ruptures remain the same but the Segmented
ruptures now have a GEA12 branch (same as 2014/2018) and a  GEA17 branch. Sections have
also been split into 50-km subsections (e.g., 2014/2018 Section 1 is now 7 50-km sections
[1-1:1-7], Section 2 is now 3 50-km sections [2-1:2-3], Section 3 is now 3 50-km sections
[3-1:3-3], and Section 4 is now 9 50-km sections [4-1:4-9]). All section IDs have been
updated, and new rupture IDs have been assigned.

| Sections                            | 1-1:1-7    | 2-1:2-3    | 3-1:3-3    | 4-1:4-9    |                            |
|:----------------------------------- |:----------:|:----------:|:----------:|:----------:|:--------------------------:|
| Bottom                              | 3200..3206 | 3207..3209 | 3210..3212 | 3213..3221 |                            |
| Middle                              | 3230..3236 | 3237..3239 | 3240..3242 | 3243..3251 |                            |
| Top                                 | 3260..3266 | 3267..3269 | 3270..3272 | 3273..3281 |                            |
| __Ruptures__                        |            |            |            |            | [ Bottom, Middle, Top ]    |
| Full, GEA12-A, 1-1:4-9              | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◆        | ◆◆◆◆◆◆◆◆◆  | [ 100002, 100001, 100000 ] |
| Segmented, GEA12, D, 1-1:1-7        | ◆◆◆◆◆◆◆    |            |            |            | [ 100302, 100301, 100300 ] |
| Segmented, GEA12, C, 1-1:2-3        | ◆◆◆◆◆◆◆    | ◆◆◆        |            |            | [ 100202, 100201, 100200 ] |
| Segmented, GEA12, B, 1-1:3-3        | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◆        |            | [ 100102, 100101, 100100 ] |
| Segmented, GEA12, Northern, 4-1:4-9 |            |            |            | ◆◆◆◆◆◆◆◆◆  | [ 100402, 100401, 100400 ] |
| Segmented, GEA17, E, 1-3:1-4        | ◇◇◆◆◇◇◇    |            |            |            | [ 100802, 100801, 100800 ] |
| Segmented, GEA17, D, 1-1:1:6        | ◆◆◆◆◆◆◇    |            |            |            | [ 100702, 100701, 100700 ] |
| Segmented, GEA17, C', 1-1:2-3       | ◆◆◆◆◆◆◆    | ◆◆◆        |            |            | [ 100202, 100201, 100200 ] |
| Segmented, GEA17, C, 1-1:3-2        | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◇        |            | [ 100602, 100601, 100600 ] |
| Segmented, GEA17, B, 1-1:4-3        | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◆        | ◆◆◆◇◇◇◇◇◇  | [ 100502, 100501, 100500 ] |
| Segmented, GEA17, F, 3-3:4-4        |            |            | ◇◇◆        | ◆◆◆◆◇◇◇◇◇  | [ 100902, 100901, 100900 ] |
| Unsegmented, GEA12-B, 1-1:3-3       | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◆        |            | [ 100102, 100101, 100100 ] |
| Unsegmented, GEA12-A, 1-1:4-9       | ◆◆◆◆◆◆◆    | ◆◆◆        | ◆◆◆        | ◆◆◆◆◆◆◆◆◆  | [ 100002, 100001, 100000 ] |

| Sections             | 1    | 2    | 3    | 4    |                         |
|:-------------------- |:----:|:----:|:----:|:----:|:-----------------------:|
| Bottom               | 3200..3206 | 3207..3209 | 3210..3212 | 3213..3221 | |
| Middle               | 3230..3236 | 3237..3239 | 3240..3242 | 3243..3251 | |
| Top                  | 3260..3266 | 3267..3269 | 3270..3272 | 3273..3281 | |
| __Ruptures__         |      |      |      |      | [ Bottom, Middle, Top ] |
| Full, 1-2-3-4        | ◆    | ◆    | ◆    | ◆    | [ 3170, 3171, 3172 ]    |
| Segmented, 1         | ◆    |      |      |      | [ 3110, 3111, 3112 ]    |
| Segmented, 1-2       | ◆    | ◆    |      |      | [ 3150, 3151, 3152 ]    |
| Segmented, 1-2-3     | ◆    | ◆    | ◆    |      | [ 3160, 3161, 3162 ]    |
| Segmented, 4         |      |      |      | ◆    | [ 3140, 3141, 3142 ]    |
| Unsegmented, 1-2-3   | ◆    | ◆    | ◆    |      | [ 3160, 3161, 3162 ]    |
| Unsegmented, 1-2-3-4 | ◆    | ◆    | ◆    | ◆    | [ 3170, 3171, 3172 ]    |

## Cluster Model (2023)

Table of clusters and relevant section ranges. The various combinations of ruptures below, each
composed of a minimum of 4 subsections, yields exactly one cluster section boundary between
sections 1-4 and 4-6.

Note on cluster- and rupture-IDs as implemented for bottom, middle, and top cluster
model files - within a `cluster-set.json` file, the cluster and rupture IDs are
incremental. Bottom IDs are in the range 3300 to 3330, middle IDs are in the range 3340 to 3370,
and top IDs are in the range 3380 to 3410.

| Cluster | 2     | 3a    | 3b    | 4a      | 4b      | 5         | 7a            | 7b            | 8               |
|:-------:|:-----:|:-----:|:-----:|:-------:|:-------:|:---------:|:-------------:|:-------------:|:---------------:|
| Section | 11-11 | 8-6-8 | 7-8-7 | 6-6-5-5 | 5-5-6-6 | 4-5-4-5-4 | 3-3-3-3-3-3-4 | 4-3-3-3-3-3-3 | 2-3-3-3-3-3-3-2 |
| 4-9     | ◇     | ◆     | ◆     | ◇       | ◇       | ◆         | ◆             | ◆             | ◆               |
| 4-8     | ◇     | ◆     | ◆     | ◇       | ◇       | ◆         | ◆             | ◆             | ◆               |
| 4-7     | ◇     | ◆     | ◆     | ◇       | ◇       | ◆         | ◆             | ◆             | ◇               |
| 4-6     | ◇     | ◆     | ◆     | ◇       | ◇       | ◆         | ◆             | ◇             | ◇               |
| 4-5     | ◇     | ◆     | ◆     | ◇       | ◇       | ◇         | ◇             | ◇             | ◇               |
| 4-4     | ◇     | ◆     | ◆     | ◆       | ◇       | ◇         | ◇             | ◇             | ◆               |
| 4-3     | ◇     | ◆     | ◆     | ◆       | ◆       | ◇         | ◇             | ◆             | ◆               |
| 4-2     | ◇     | ◆     | ◇     | ◆       | ◆       | ◇         | ◆             | ◆             | ◆               |
| 4-1     | ◇     | ◇     | ◇     | ◆       | ◆       | ◇         | ◆             | ◆             | ◇               |
| 3-3     | ◇     | ◇     | ◇     | ◆       | ◆       | ◆         | ◆             | ◇             | ◇               |
| 3-2     | ◇     | ◇     | ◇     | ◇       | ◆       | ◆         | ◇             | ◇             | ◇               |
| 3-1     | ◆     | ◇     | ◇     | ◇       | ◆       | ◆         | ◇             | ◇             | ◆               |
| 2-3     | ◆     | ◇     | ◇     | ◇       | ◇       | ◆         | ◇             | ◆             | ◆               |
| 2-2     | ◆     | ◇     | ◇     | ◇       | ◇       | ◇         | ◆             | ◆             | ◆               |
| 2-1     | ◆     | ◆     | ◇     | ◇       | ◇       | ◇         | ◆             | ◆             | ◇               |
| 1-7     | ◆     | ◆     | ◆     | ◇       | ◇       | ◇         | ◆             | ◇             | ◇               |
| 1-6     | ◆     | ◆     | ◆     | ◆       | ◇       | ◇         | ◇             | ◇             | ◇               |
| 1-5     | ◆     | ◆     | ◆     | ◆       | ◆       | ◇         | ◇             | ◇             | ◆               |
| 1-4     | ◆     | ◆     | ◆     | ◆       | ◆       | ◆         | ◇             | ◆             | ◆               |
| 1-3     | ◆     | ◆     | ◆     | ◆       | ◆       | ◆         | ◆             | ◆             | ◆               |
| 1-2     | ◆     | ◆     | ◆     | ◆       | ◆       | ◆         | ◆             | ◆             | ◇               |
| 1-1     | ◆     | ◆     | ◆     | ◆       | ◆       | ◆         | ◆             | ◆             | ◇               |

## Rupture Reference (2014, 2018)

The table below lists the IDs used for rupture sets on different branches of the 2014 & 2018
Cascadia logic tree (based on the GEA12 model).

| Sections                | 1    | 2    | 3    | 4    |                         |
|:----------------------- |:----:|:----:|:----:|:----:|:-----------------------:|
| Bottom                  | 3110 | 3120 | 3130 | 3140 |                         |
| Middle                  | 3111 | 3121 | 3131 | 3141 |                         |
| Top                     | 3112 | 3122 | 3132 | 3142 |                         |
| __Ruptures__            |      |      |      |      | [ Bottom, Middle, Top ] |
| Full, A, 1-2-3-4        | ◆    | ◆    | ◆    | ◆    | [ 3170, 3171, 3172 ]    |
| Segmented, D, 1         | ◆    |      |      |      | [ 3110, 3111, 3112 ]    |
| Segmented, C, 1-2       | ◆    | ◆    |      |      | [ 3150, 3151, 3152 ]    |
| Segmented, B, 1-2-3     | ◆    | ◆    | ◆    |      | [ 3160, 3161, 3162 ]    |
| Segmented, Northern, 4  |      |      |      | ◆    | [ 3140, 3141, 3142 ]    |
| Unsegmented, B, 1-2-3   | ◆    | ◆    | ◆    |      | [ 3160, 3161, 3162 ]    |
| Unsegmented, A, 1-2-3-4 | ◆    | ◆    | ◆    | ◆    | [ 3170, 3171, 3172 ]    |

## Logic Tree Weight Adjustments

__bot:__ 1.856 instead of 1.8534  
__mid:__ 1.8528 instead of 1.8534  
__top:__ 1.856 instead of 1.8534  
