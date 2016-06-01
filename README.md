<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Linux Build Status](https://travis-ci.org/SMAC-Group/gmwm.svg?branch=master)](https://travis-ci.org/SMAC-Group/gmwm)[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/gmwm)](http://www.r-pkg.org/pkg/gmwm)[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/gmwm)](http://cran.r-project.org/package=gmwm)

`gmwm` R Package
================

This repository holds the Generalized Method of Wavelet Moments (GMWM) R package. This estimation technique uses the wavelet variance in a moment-matching spirit to estimate parameters of time series models such as ARMA or state-space models. A robust version of the GMWM is also implemented in this package together with a tailored made method for inertial sensors calibration, which typically deals with very large sample sizes.

Below are examples of the capabilities of the `gmwm` package.

To start, let's generate some data:

``` r
## Data generation ##
# Specify model
m = AR1(phi = .98, sigma2 = .01) + WN(sigma2 = 1)

# Generate Data
d = gen.gts(m, 10000)
```

Once we have data, we can see what the wavelet variance looks like for the data with the classical and robust wavelet variances.

``` r
# Calculate the classical wavelet variance with the Haar filter
wv.classical = wvar(d)

# Plot the data
plot(wv.classical)

# Calculate robust wavelet variance
wv.robust = wvar(d, robust = TRUE, eff = 0.6)

# Compare both versions
compare.wvar(wv.classical, wv.robust)
```

Now, let's try to estimate it with specific (e.g. user supplied) and guessed (e.g. program generated) parameters.

``` r
## Estimation Modes ##

# Use a specific initial starting value
o.specific = gmwm.imu(AR1(phi=.98,sigma2=.05) + WN(sigma2=.95), data = d)

# Let the program guess a good starting value
o.guess = gmwm.imu(AR1()+WN(), data = d)
```

To run inference or view the parameter estimates, we do:

``` r
## View Model Info ##

# Standard summary
summary(o.specific)

# View with asymptotic inference
summary(o.specific, inference = T)

# View with bootstrapped inference
summary(o.specific, inference = T, bs.gof = T)
```

Alternatively, we can let the program try to figure out the best model for the data using the Wavelet Information Criteria (WIC):

``` r
## Model selection ##

# Separate Models - Compares 2*AR1() and AR1() + WN() under common model 2*AR1() + WN()
# Note: This function created a shared model (e.g. 2*AR1() + WN()) if not supplied to obtain the WIC. 
ms.sep = rank.models(AR1()+WN(), 2*AR1(), data = d, model.type="imu")

# Nested version - Compares AR1() + WN(), AR1(), WN()
ms.nested = rank.models(AR1()+WN(), data = d, nested = TRUE, model.type = "imu")

# Bootstrapped Optimism
ms.bs = rank.models(AR1()+WN(), WN(), data = d, bootstrap = TRUE, model.type = "imu")

# See automatic selection fit
plot(ms.sep)

# View model picked:
summary(ms.sep)
```

Last, but certainly not least, we can also approximate a contaminated sample with robust methodology:

``` r
## Data generation ##
# Specify model
model = AR1(phi = .99, sigma2 = .01) + WN(sigma2 = 1)

# Generate Data
set.seed(213)
N = 1e3
sim.ts = gen.gts(model, N)

# Contaminate Data
cont.eps = 0.01
cont.num = sample(1:N,round(N*cont.eps))
sim.ts[cont.num,] = sim.ts[cont.num,] + rnorm(round(N*cont.eps),0,sqrt(100))

# Plot the data
plot(sim.ts)

# Classical Wavelet Variance
wv.classic = wvar(sim.ts)

# Robust Wavelet Variance
wv.robust = wvar(sim.ts, robust = TRUE, eff = 0.6)

# Plot the Classical vs. Robust WV
compare.wvar(wv.classic, wv.robust, split = FALSE)

# Run robust estimation
o = gmwm.imu(model, sim.ts, robust = TRUE, eff = 0.6)

# Robust information
summary(o)
```

Install Instructions
====================

To install the `gmwm` package, there are three options: CRAN (Stable), GitHub (Developmental), or SMAC (stable - offline).

Recommended R Interface
-----------------------

We firmly recommend that any users of this package use the [RStudio IDE](https://www.rstudio.com/products/rstudio/download/) over the default R GUI.

Installing the package through CRAN (Stable)
--------------------------------------------

The installation process with CRAN is the simplest

``` r
install.packages("gmwm")
```

Installing the package this way gives you access to stable features. Furthermore, the installation itself does not require a compiler or preinstalling any dependencies. However, we are limited to updating the package on CRAN to once every month. Thus, there may be some lag between when features are developed and when they are available on this version.

Installing the package through GitHub (Developmental)
-----------------------------------------------------

For users who are interested in having the latest and greatest developments withing wavelets or GMWM methodology, this option is ideal. Though, there is considerably more work that a user must do to have a stable version of the package. **The setup to obtain the development version is platform dependent.**

Specifically, one **must** have a compiler installed on your system that is compatible with R.

For help on obtaining a compiler consult:

-   [OS X](http://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-os-x/)
-   [Windows](https://cran.r-project.org/bin/windows/Rtools/)

Depending on your operating system, further requirements exist such as:

**OS X**

Some user report the need to use X11 to suppress shared library errors. To install X11, visit [xquartz.org](http://www.xquartz.org/)

**Linux**

Both curl and libxml are required.

For **Debian** systems, enter the following in terminal:

``` bash
sudo apt-get install curl libcurl3 libcurl3-dev libxml2 libxml2-dev
```

For **RHEL** systems, enter the following in terminal:

``` bash
sudo yum install curl curl-devel libxml2 libxml2-dev
```

**All Systems**

With the system dependency taken care of, we continue on by installing the R specific package dependencies and finally the package itself by doing the following in an R session:

``` r
# Install dependencies
install.packages(c("RcppArmadillo","ggplot2","reshape2","devtools","knitr","rmarkdown"))

# Install the package from github
devtools::install_github("SMAC-Group/gmwm")
```

Installing the package from SMAC-Group.com (Stable - offline)
-------------------------------------------------------------

Lastly, we will be offering a source .tar that is able to be install offline - after being downloaded - on the [smac-group.com](http://www.smac-group.com) website.

``` r
# Install the dependencies
install.packages(c("RcppArmadillo","ggplot2","scales","devtools","knitr","rmarkdown"))

# Local installation
setwd("path_to_file_GMWM_2.0.0.tar.gz")
install.packages("GMWM", repos = NULL, type="source")
```

Supplementary data package
--------------------------

To test the package performance on real-world data that is *stationary* or work with some of the examples, you will need to download and install the `imudata` R package.

To do so, please use the following installation method within the `gmwm` R package:

``` r
gmwm::install_imudata()
```

For more information about the `imudata` package, see the [repository](https://github.com/SMAC-Group/imudata).

Licensing
=========

The license this source code is released under is the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA). In some cases, the GPL license does apply. However, in the majority of the cases, the license in effect is the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA) as the computational code is heavily dependent on armadilllo, which has an MIT license that enables us to recast our code to the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA). See the LICENSE file for full text. Otherwise, please consult [TLDR Legal](https://tldrlegal.com/license/creative-commons-attribution-noncommercial-sharealike-(cc-nc-sa)) or [CC](https://creativecommons.org/licenses/by-nc-sa/4.0/#) which will provide a synopsis of the restrictions placed upon the data and code. Please note, this does NOT excuse you from talking about licensing with a lawyer!
