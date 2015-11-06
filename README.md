# `gmwm` R Package
This repository holds the Generalized Method of Wavelet Moments (GMWM) R package. This estimation technique uses the wavelet variance in a moment-matching spirit to estimate parameters of time series models such as ARMA or state-space models. A robust version of the GMWM is also implemented in this package together with a tailored made method for inertial sensors calibration, which typically deals with very large sample sizes. 

Below are examples of the capabilities of the `gmwm` package.

To start, let's generate some data:
```r
## Data generation ##
# Specify model
m = AR1(phi=.99,sigma2=.01) + WN(sigma2=1)

# Generate Data
d = gen.gts(m, 10000)
```

Once we have data, we can see what the wavelet variance looks like for the data with the classical and robust wavelet variances.

```r
# Calculate the classical wavelet variance with the Haar filter
wv.classical = wvar(d)

# Plot the data
plot(wv.classical)

# Calculate robust wavelet variance
wv.robust = wvar(d, robust = T, eff = 0.6)

# Compare both versions
compare.wvar(wv.classical, wv.robust)
```

Now, let's try to estimate it with specific (e.g. user supplied) and guessed (e.g. program generated) parameters.

```r
## Estimation Modes ##

# Use a specific initial starting value
o.specific = gmwm.imu(AR1(phi=.95,sigma2=.05) + WN(sigma2=.95), d)

# Let the program guess a good starting value
o.guess = gmwm.imu(AR1()+WN(), d)
```

To run inference or view the parameter estimates, we do:
```r
## View Model Info ##

# Standard summary
summary(o.specific)

# View with inference
summary(o.specific, inference = T)

# Add bootstrapping
summary(o.specific, inference = T, bs.gof = T)
```

Alternatively, we can let the program try to figure out the best model for the data:
```r
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

Last, but certainly not least, we can also approximate a contaminated sample with robust methodology:
```r
## Data generation ##
# Specify model
model = AR1(phi=.99,sigma2=.01) + WN(sigma2=1)

# Generate Data
data = gen.gts(model, N)

# Contaminate Data
cont.alpha = 0.01
cont.num = round(N*cont.alpha)
data$data[sample(1:N,cont.num),] = rnorm(cont.num)

# Robust wavelet variance
wv.robust = wvar(data, robust = T, eff = 0.6)

# Plot the robust WV
plot(wv.robust)

# Run robust estimation
o = gmwm.imu(model, data, robust = T, eff = 0.6)

# Robust information
summary(o)
```


# Install Instructions (All platforms)
To install the `gmwm` package, there are three options: CRAN (stable), GitHub (Developmental), or SMAC (stable - offline).

The installation process with CRAN is the simplest
```r
install.packages("gmwm")
```

Prior to installing with `devtools`, please make sure to have a compiler installed on your system that is compatible with R.

If you have a compiler already installed, then continue on by installing the package dependencies and finally the package itself by doing the following: 

```r
# Install dependencies
install.packages(c("RcppArmadillo","ggplot2","gridExtra","reshape2","devtools"))

# Install the package from github
devtools::install_github("SMAC-Group/gmwm")
```

Lastly, we will be offering a source .tar that is able to be install offline - after being downloaded - on the [smac-group.com](http://www.smac-group.com) website.

```r
# Install the dependencies
install.packages(c("RcppArmadillo","ggplot2", "scales", "gridExtra","devtools"))

# Local installation
setwd("path_to_file_GMWM_0.13.0.tar.gz")
install.packages("GMWM", repos = NULL, type="source")
```

## Supplementary data package

To test the package performance on real-world data that is *stationary* or work with some of the examples, you will need to download and install the `imudata` R package.

To do so, please use the following installation method within the `gmwm` R package:

```r
gmwm::install_imudata()
```

For more information about the `imudata` package, see the [repository](https://github.com/SMAC-Group/imudata).

# Licensing
The license this source code is released under is the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA). In some cases, the GPL license does apply. However, in the majority of the cases, the license in effect is the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA) as the computational code is heavily dependent on armadilllo, which has an MIT license that enables us to recast our code to the Creative Commons Attribution NonCommercial ShareAlike (CC-NC-SA). See the LICENSE file for full text. Otherwise, please consult [TLDR Legal](https://tldrlegal.com/license/creative-commons-attribution-noncommercial-sharealike-(cc-nc-sa)) or [CC](https://creativecommons.org/licenses/by-nc-sa/4.0/#) which will provide a synopsis of the restrictions placed upon the data and code. Please note, this does NOT excuse you from talking about licensing with a lawyer!
