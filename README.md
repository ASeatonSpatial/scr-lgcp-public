# Spatial capture-recapture with penalized regression splines

This repository contains code and data to accompany the following paper:

* Seaton, A. E., Borchers, D. L., Groenenberg, M., and Stevenson, B. C. (in submission) Spatial capture-recapture with penalized regression splines to flexibly model wildlife density and distribution.

## Dependencies

The following R packages are required to run the code in this repository, and are available on CRAN: `fields`, `mgcv`, `Rcpp`, `RColorBrewer`, `sp`, `TMB`, and `viridis`.

The following R packages are also required, but are not available on CRAN:
* `acre`: Install by following the instructions [here](https://github.com/b-steve/ascr).
* `INLA`: Install by following the instructions [here](https://www.r-inla.org/download/index.html).

Note that both `Rcpp` and `TMB` also require a C++ compiler. Windows users can obtain one by installing [`RTools`](https://cran.r-project.org/bin/windows/Rtools/).

