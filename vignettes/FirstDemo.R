## @knitr library_include
library(GMWM)

## @knitr gen_process
# Set seed for reproducibility
set.seed(1)

# Length of the time series
n = 100000

# Simulate AR(1) + WN
xt = gen.gts(AR1(phi=.99, sigma2 = 0.01) + WN(sigma2=1),n)

## @knitr wv
wv = wvar(xt)
plot(wv)

## @knitr modelTS
TS.mod = AR1() + WN()

## @knitr GMWM
model = gmwm(TS.mod, data = xt)

## @knitr results
results = matrix(c(round(model$estimate,4), c(0.99,0.01,1)), 3, 2)
dimnames(results)[[1]] = c("phi","sig2AR1","sig2WN")
dimnames(results)[[2]] = c("Estim.","True")
results

## @knitr GMWMplot
plot(model)

## @knitr imuData
data(imu)
head(imu)

## @knitr imuPlot separate
wv.imu = imu2WV(imu)
plot(wv.imu)

## @knitr imuPlot combined
plot(wv.imu, separate=FALSE)

## @knitr imuModel
# Define model
TS.mod.imu = 3*AR1()

# Compute GMWM estimator
model.imu = gmwm(TS.mod.imu, data = imu[,1])

## @knitr imuModel summary
summary(model.imu)
plot(model.imu)

## @knitr newdata
# Add more...