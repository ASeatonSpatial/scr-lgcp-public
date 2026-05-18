# Spatial capture-recapture with penalized regression splines

This repository contains code and data to accompany the following paper:

* Seaton, A. E., Borchers, D. L., Groenenberg, M., and Stevenson, B. C. (in submission) Spatial capture-recapture with penalized regression splines to flexibly model wildlife density and distribution.

## Dependencies

The following R packages are required to run the code in this repository, and are available on CRAN: `fields`, `mgcv`, `Rcpp`, `RColorBrewer`, `sp`, `TMB`, and `viridis`.

The following R packages are also required, but are not available on CRAN:
* `acre`: Install by following the instructions [here](https://github.com/b-steve/ascr).
* `INLA`: Install by following the instructions [here](https://www.r-inla.org/download/index.html).

Note that both `Rcpp` and `TMB` also require a C++ compiler. Windows users can obtain one by installing [`RTools`](https://cran.r-project.org/bin/windows/Rtools/).

## Gibbon case study

Code for the gibbon case study is available in `gibbons.R`. The data are available in `gibbon-data.RData`.

## Black bear case study

We do not yet have permission to publicly share the black bear data. We are hoping to obtain permission soon, at which point we will upload the data along with our code.

## Simulation study

Code to set up our simulation study, then fit all models to one data set (i.e., a single simulation iteration) can be found in `sim.R`. We carried out our simulation study using parallel processing on high-performance computing facilities provided by New Zealand's eScience Infrastructure.
