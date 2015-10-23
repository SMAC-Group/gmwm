# `gmwm` R Package
This repository holds the Generalized Method of Wavelet Moments (GMWM) Statistical Methods R package. The methods within this package are able to be used to estimate parameters of latent time series models within two different disciplines: Inertial Measurement Units for Engineers and State-Space Models for Statisticians.

Example workflow:
```r
## Data generation ##
# Specify model
m = AR1(phi=.99,sigma2=.01) + WN(sigma2=1)

# Generate Data
d = gen.gts(m, 10000)

## Estimation Modes ##

# Use a specific initial starting value
o.specific = gmwm.imu(AR1(phi=.95,sigma2=.05) + WN(sigma2=.95), d)

# Let the program guess a good starting value
o.guess = gmwm.imu(AR1()+WN(), d)

## View Model Info ##

# Standard summary
summary(o.specific)

# View with inference
summary(o.specific, inference = T)

# Add bootstrapping
summary(o.specific, inference = T, bs.gof = T)

## Model selection ##

# Separate Models - Compares 2*AR1() + WN(), 2*AR1(), AR1() + WN()
ms.sep = rank.models(d,models=list(AR1()+WN(),2*AR1()), model.type="imu")

# Nested version - Compares AR1() + WN(), AR1(), WN()
ms.nested = rank.models(d,models=list(AR1()+WN()), nested = T, model.type="imu")

# Bootstrapped Optimism
ms.bs = rank.models(d,models=list(AR1()+WN(),WN()), bs.optimism = T, model.type="imu")

# See automatic selection fit
plot(ms.sep)

# View model picked:
summary(ms.sep)
```

# Install Instructions (All platforms)
To install the `gmwm` package, there are two options: CRAN or github.

Prior to installing with devtools, please make sure to have a compiler installed on your system that is compatible with R.

If you have a compiler already installed, then continue on by installing the package dependencies and finally the package itself by doing the following: 

```r
# Install dependencies
install.packages(c("RcppArmadillo","ggplot2","gridExtra","reshape2","devtools"))

# Install the package from github
devtools::install_github("SMAC-Group/imudata")
```

## Supplementary data package

To test the package performance on real-world data that is *stationary* or work with some of the examples, you will need to download and install the `imudata` R package.

To do so, please use the following installation method within the `gmwm` R package:

```r
gmwm::install_imudata()
```

Lastly, we will be offering a source .tar that is able to be install offline - after being downloaded - on the [smac-group.com](http://www.smac-group.com) website.

```r
# Install the dependencies
install.packages(c("RcppArmadillo","ggplot2", "scales", "gridExtra","devtools"))

# Local installation
setwd(“path_to_file_GMWM_0.13.0.tar.gz”)
install.packages(“GMWM", repos = NULL, type="source")
```

# Licensing
The license this data is released under is the Q public license. See the LICENSE file for full text. Otherwise, please consult [TLDR Legal](https://tldrlegal.com/license/q-public-license-1.0-%28qpl-1.0%29) which will provide a synopsis of the restrictions placed upon the data and code.
