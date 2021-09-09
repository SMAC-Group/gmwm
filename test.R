library(gmwm)
n = 2.5*10^5
model = WN(sigma2 = 1) + AR(phi = 0.999, sigma2 = 10^(-5))
Xt = gen_gts(n = n, model = model)

# Standard GMWM using WV as input
wv_Xt = wvar(Xt)
plot(wv_Xt)
estim = gmwm(model, wv_Xt)
estim
plot(estim)

# Standard GMWM using data as input
estim = gmwm(model, Xt)
estim
plot(estim)

# Robust GMWM using WV as input
wv_Xt = wvar(Xt, robust = TRUE)
plot(wv_Xt)
estim = gmwm_wvar(model, wv_Xt)
estim
plot(estim)

# Robust GMWM using data as input
estim = gmwm(model, Xt, robust = TRUE)
estim
plot(estim)

# Testing removing some scales
wv_Xt = wvar(Xt)
plot(wv_Xt)

# Making the first scales looks weird!
wv_Xt$variance[1:2] = wv_Xt$variance[1:2] + c(5, 1)
wv_Xt$ci_low[1:2] = wv_Xt$ci_low[1:2] + c(5, 1)
wv_Xt$ci_high[1:2] = wv_Xt$ci_high[1:2] + c(5, 1)
plot(wv_Xt)

estim = gmwm(model, wv_Xt)
estim
plot(estim)


estim = gmwm(model, wv_Xt, remove_scales = 1:2)
estim
plot(estim)

# Use custom Omega matrix
wv_Xt = wvar(Xt)
estim = gmwm(model, wv_Xt)
estim
plot(estim)

estim = gmwm(model, wv_Xt, Omega = diag(rep(1, wv_Xt$J)))
estim
plot(estim)


