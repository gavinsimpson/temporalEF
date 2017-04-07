# temporalEF

#### Released version
[![CRAN version](http://www.r-pkg.org/badges/version/temporalEF)](http://cran.rstudio.com/web/packages/temporalEF/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/temporalEF)](http://cran.rstudio.com/web/packages/temporalEF/index.html)

#### Build status
[![Build Status](https://travis-ci.org/gavinsimpson/temporalEF.svg?branch=master)](https://travis-ci.org/gavinsimpson/temporalEF)  [![Build status](https://ci.appveyor.com/api/projects/status/hc8dbxrim2nj3c1i/branch/master)](https://ci.appveyor.com/project/gavinsimpson/temporalEF/branch/master)  [![codecov](https://codecov.io/gh/gavinsimpson/temporalEF/branch/master/graph/badge.svg)](https://codecov.io/gh/gavinsimpson/temporalEF)

## What is temporalEF?

**temporalEF** is an R package that computes eigenfunctions for use in
explaining patterns in multivariate time series. The basic idea is taken
from recent developments in the field of spatial analysis in multivariate
data sets.

Principal coordinates of neighbour matrices (PCNM) and the more recent
Asymmetric eigenvector maps (AEM) were designed to model variation at
multiple spatial scales. **temporalEF** utilises these methods by taking
the time direction as a one-dimensional spatial process. The resulting
eigenfunctions represent sinusoidal functions of varying periodicity.

The resulting eigenfunctions can be used in constrained ordinations such
as redundancy analysis (RDA), canonical correspondence analysis (CCA),
or constrained analysis of principal coordinates (CAP) to decompose a
species data matrix into independent patterns of temporal variation or
trends.

## Licence

**temporalEF** is made available under the [GNU GPL version 2
licence](http://www.gnu.org/licenses/gpl-2.0.html).
