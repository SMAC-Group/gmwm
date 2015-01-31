## Benchmarking
library(microbenchmark)

## Waveslim package
library(waveslim)

### Support Functions ###

x = 1:5
p = c(0.1,.2,.3,.5,.6)

### Tested using vectors since C++ is set to handle vector input.

# ===Logit Inv===
microbenchmark(arma_fun=pseudo.logit.inv(x), r_fun = pseudo.logit.inv(x))

# ===Logit ===
microbenchmark(arma_fun=pseudo.logit(p), r_fun = pseudo.logit(p))


### Tested using doubles since C++ can only handle double input.
omega = 4
tau = 3

# ===Expected value DR===
microbenchmark(arma_fun=e_drift(omega, tau), r_fun = E.drift(omega, tau))

# ===Second moment DR===
microbenchmark(arma_fun=m2_drift(omega, tau), r_fun = m2.drift(omega, tau))

# ===Variance DR===
microbenchmark(arma_fun=var_drift(omega, tau), r_fun = var.drift(omega, tau))


### Generation Code ###
N = 1000
phi = .3
slope = 4
sig2 = 1.2
x.sim = 1:N

# ===AR 1 Generation===
microbenchmark(arma_fun=gen_ar1(N,phi,sig2), r_fun = doAR1(N,phi,sig2))

# ===WN Generation===
microbenchmark(arma_fun=gen_white_noise(N,sig2), r_fun = doWhiteNoise(N, sig2))

# ===RW Generation===
microbenchmark(arma_fun=gen_drift(N,slope), r_fun = doAR1(N,phi,sig2))



### x to WV ###
N = 1000
sig2 = 1.5
tau = floor(log(N,2)) #c(1,2,3,4,5,floor(log(N,2)))
phi=.37

## Note, scales has been dropped from arma output since tau is already active.

# ===WN to WV===
microbenchmark(arma_fun=wn_to_wv(sig2, tau), r_fun = White.Noise2WV(sig2, tau))

# ===RW to WV===
microbenchmark(arma_fun=rw_to_wv(sig2, tau), r_fun = Random.Walk2WV(sig2, tau))

# ===DR to WV===
microbenchmark(arma_fun=dr_to_wv(phi,sig2), r_fun = Drift2WV(sig2,tau))

# ===AR1 to WV===
microbenchmark(arma_fun=ar1_to_wv(phi,sig2,tau), r_fun = AR1.2.WV(phi,sig2,tau))

### Code ported from waveslim
N = 2^8
signal = rnorm(N)
nb.level = floor(log(N,2))
strWavelet = "haar"
haar=c(1/sqrt(2),1/sqrt(2))
p = 0.025 #for CI

### === qmf
microbenchmark(arma_fun=qmf(haar))

### === haar_filter construction
microbenchmark(arma_fun=haar_filter())

### === select_filter
microbenchmark(arma_fun=select_filter("haar"), r_fun = wave.filter("haar"))

### === dwt ===
microbenchmark(arma_fun=dwt_arma(signal,"haar",nb.level), r_fun = dwt(signal,"haar",nb.level))

### === modwt ===
microbenchmark(arma_fun=modwt_arma(signal,"haar",nb.level), r_fun = modwt(signal,"haar",nb.level))

#-- Storing Signals
signal.modwt = modwt(signal, strWavelet, nb.level)
signal.modwt.bw = brick.wall(signal.modwt, strWavelet)
signal_modwt = modwt_arma(signal, strWavelet, nb.level)
signal_modwt_bw = brick_wall(signal_modwt, haar_filter())

### === brick_wall ===
microbenchmark(arma_fun=brick_wall(signal_modwt, haar_filter()), r_fun = brick.wall(signal.modwt, strWavelet))

### === Wave variance
microbenchmark(arma_fun=wave_variance(signal_modwt_bw, "eta3", p), r_fun = wave.variance(signal.modwt.bw, p = p, type = "eta3"))

fft.in = signal.modwt[[1]]
fft_in = signal_modwt[[1]]
### === my.acf (e.g. fourier transforms) === We lose here =( Note fft and ifft are not optimized per http://arma.sourceforge.net/docs.html#fft
microbenchmark(arma_fun=dft_acf(fft_in), r_fun = my.acf(fft.in))

### === waveletVariance #Need to load demo file for functions
microbenchmark(arma_fun=wavelet_variance_arma(signal,strWavelet="haar",compute_v="diag"),
               r_fun = waveletVariance(signal, strWavelet="haar", verbose=FALSE, compute.v="diag"))

### Objective Functions


### Using the simGMWM function

### Using the GMWM function






