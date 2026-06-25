# ucla_plha

Welcome to the documentation of the **ucla_plha** Python package for performing
probabilistic liquefaction hazard analysis (PLHA). We developed it as a lightweight
and computationally efficient engine for performing PLHA calculations, in order to
quantify the impacts of sea level rise on the annual rate of liquefaction in coastal
areas of California. The work presented here was funded through the University of
California Climate Action Initiative as a Seed Award under Application ID R02CP7031.

We are making this software available because we believe it may be useful to other
researchers studying liquefaction, or any other problem involving earthquake-induced
damage where component fragility can be modeled as a lognormal distribution function.
The code can be run on a single core and has also been configured to run on HPC
resources at the Texas Advanced Computing Center ([TACC](https://tacc.utexas.edu)).
Source code is available via [GitHub](https://github.com/meerakota/ucla_plha).

```{figure} images/Figure1.png
:alt: A sketch showing ocean waves meeting the beach with a sandy berm separating the beach from infrastructure systems including pipelines, buildings, and transportation. Increases in groundwater level from sea level rise are illustrated by black dotted lines. Liquefaction is more likely to occur and damage infrastructure systems when groundwater rises.
:width: 800px
:align: center

Sea level rise will increase groundwater elevation, thereby increasing the frequency
of earthquake-induced liquefaction.
```

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 🚀 Get started
:link: get_started
:link-type: doc
Install the package and run your first hazard analysis from a JSON config.
:::

:::{grid-item-card} 📖 Functions
:link: functions
:link-type: doc
Reference for every function in the `ucla_plha` package.
:::

:::{grid-item-card} 🧪 Examples
:link: examples
:link-type: doc
Worked config files and a Jupyter notebook covering the full workflow.
:::

:::{grid-item-card} ✅ Validation
:link: validation
:link-type: doc
PSHA, ground motion, liquefaction, and disaggregation checks against published benchmarks.
:::

::::

```{toctree}
:maxdepth: 2
:caption: Documentation
:hidden:

get_started
overview
config_format
examples
validation
functions
```
