# Define functions used in gmwm fct when processes can be expressed linearly w/ parameters but that do not need to be called by the user

# # test equation of the theoretical wavelet variance for stochastic processes
# # that can be expressed linearly with respect to the parameters
# 
# # load libraries
# library(avar)
# library(wv)
# library(simts)
# 
# # WN
# true_sigma2 = 10
# Xt = rnorm(1000, sd = sqrt(true_sigma2))
# my_wvar = wvar(Xt)
# plot(my_wvar)
# lines(x = my_wvar$scales, y = true_sigma2/my_wvar$scales, col = "orange", lwd= 2, type ="b")
# 
# # RW
# true_gamma2 = 10
# my_mod =  simts::RW(gamma2 = true_gamma2)
# Xt = gen_gts(n = 1000, model = my_mod)
# my_wvar = wvar(Xt)
# plot(my_wvar)
# lines(x = my_wvar$scales, y = ( (my_wvar$scales^2+2) * true_gamma2) / (12*my_wvar$scales) , col = "orange", lwd= 2, type ="b")
# 
# # QN
# true_q2 = 10
# my_mod =  simts::QN(q2 = true_q2)
# Xt = gen_gts(n = 1000, model = my_mod)
# my_wvar = wvar(Xt)
# plot(my_wvar)
# lines(x = my_wvar$scales, y = (6 * true_q2) / my_wvar$scales^2 , col = "orange", lwd= 2, type ="b")
# 
# # DR
# true_omega = 10
# my_mod =  simts::DR(omega = true_omega)
# Xt = gen_gts(n = 1000, model = my_mod)
# my_wvar = wvar(Xt)
# plot(my_wvar)
# lines(x = my_wvar$scales, y = (my_wvar$scales^2*true_omega^2)/16 , col = "orange", lwd= 2, type ="b")


# Define functions used in gmwm fct when processes can be expressed linearly w/ parameters but that do not need to be called by the user

# ## TO REMOVE???
# return_Omega = function(y){
#   wv_y = gmwm::wvar(y)
#   return(diag(1/(wv_y$ci_high - wv_y$ci_low)^2))
# }
# 
# ## TO REMOVE???
# return_matrix = function(model, y){
#   #define tau
#   n = length(y)
#   J = floor(log2(n))
#   J_vec = seq(J)
#   def_scales = 2^J_vec
# 
#   #define matrix X for linear in parameter processes
#   qn_p = 6 / (def_scales^2)
#   wn_p = 1/def_scales
#   rw_p = (def_scales^2+2) / (12*def_scales)
#   dr_p = (def_scales^2)/16  #note that this is linear with omega^2, not omega
#   complete_X_mat = cbind(qn_p, wn_p, rw_p, dr_p)
#   colnames(complete_X_mat) = c('QN', 'WN', 'RW', 'DR')
# 
#   #define X_mat
#   X_mat = complete_X_mat[, model$process.desc]
# 
#   #return matrix
#   return(X_mat)
# }


# TO DOCUMENT
#' @export
gmwm_wvar = function(model, wv, model.type="imu", compute.v="auto", remove_scales = NULL, Omega = NULL,
                     alpha = 0.05, seed = 1337, G = NULL, K = 1, H = 100,
                     freq = 1){
  
  # ADD SOME CHECKS HERE!
  # Define Omega 
  if (is.null(Omega)){
    Omega = diag(1/(wv$ci_high - wv$ci_low)^2)
  }
  
  # ADD SOME CHECKS HERE!
  # Remove scales that shouldn't be considered
  if (!is.null(remove_scales)){
    nb_to_remove = length(remove_scales)
    min_omega = min(diag(Omega))/10^4
    for (k in 1:nb_to_remove){
      Omega[remove_scales[k], remove_scales[k]] = min_omega
    }
  }
  
  # Do we have a valid model?
  if(!is.ts.model(model)){
    stop("`model` must be created from a `ts.model` object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  np = model$plength

  N = wv$N
  
  starting = model$starting
  
  # Input guessing
  #G=0 #modified because caused error : Error in round(x) : non-numeric argument to mathematical function"
  #if starting, g = 0
  #else is.null(G) 
  #virer le starting et mettre dans eslse
  if((is.null(G)) || !is.wholenumber(G)){
    if(N > 10000){
      G = 1e6
    }else{
      G = 20000
    }
  }
  
  # For reproducibility
  set.seed(seed)
  
  num.models = count_models(desc)
  
  # Identifiability issues
  if(any(num.models[c("DR","QN","RW","WN")]  > 1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(num.models["GM"]> 0 & num.models["AR1"] > 0){
    stop("Please either use `GM()` or `AR1()` model components. Do not mix them.")
  }
  
  # Model type issues
  model.type = tolower(model.type)
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either `ssm` or `imu`!")
  }
  
  # Verify Scales and Parameter Space
  nlevels =  wv$J
  
  if(np > nlevels){
    stop("Please supply a longer signal / time series in order to use the GMWM.",
         "This is because we need at least the same number of scales as",
         "parameters to estimate.")
  }
  
  if(wv$robust){
    np = np+1
    if(np > nlevels){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we at least need the same number of scales as parameters to estimate.")
    }
    
    if(wv$eff > 0.99){
      stop("The efficiency specified is too close to the classical case. Use `robust = FALSE`")
    }
  }
  
  # Compute fast covariance if large sample, otherwise, bootstrap.
  if(compute.v == "auto" || ( compute.v != "fast" && compute.v != "diag" &&
                              compute.v != "full" && compute.v != "bootstrap" )){
    compute.v = "fast"
  }
  
  theta = model$theta
  
  detected_gm = any(model$desc == "GM")
  
  if(detected_gm && freq == 1){
    warning("'freq' is set to 1 by default this impacts the `GM()` parameters. See ?GM for more details.")
  }
  
  # Convert from GM to AR1
  if(!starting && detected_gm){
    theta = conv.gm.to.ar1(theta, model$process.desc, freq)
  }
  
  #if sub processes are linear, fit using weighted least squares, otherwise call c++ code
  # if(all(model$process.desc %in% c('QN', 'WN', 'RW', 'DR'))){
  #   X_mat  = return_matrix(model = model, y = data)
  #   Omega  = return_Omega(data)
  #   nu_hat = gmwm::wvar(data)$variance
  #   theta_hat = solve(t(X_mat) %*% Omega %*% X_mat) %*% t(X_mat) %*% Omega %*% nu_hat
  #   colnames(theta_hat) = 'Estimates'
  #   rownames(theta_hat) = model$process.desc
  #   ci_h = gmwm::wvar(data)$ci_high
  #   ci_l = gmwm::wvar(data)$ci_low
  #   obj_teo = wv::decomp_theoretical_wv(theta = theta_hat, desc = model$process.desc, objdesc = model$obj.desc, tau = nu_hat$scales)
  #   sum_theo = if(is.vector(X_mat)){sum_theo = X_mat}else if(is.matrix(X_mat)){sum_theo = rowSums(X_mat)}
  #   out = structure(list('estimate' = theta_hat,
  #                  'init.guess' = NA,
  #                  'wv.empir' = nu_hat,
  #                  'ci.low' = ci_l,
  #                  'ci.high' = ci_h,
  #                  'orgV' = NA,
  #                  'V' = NA,
  #                  'omega' = NA,
  #                  'obj.fun' = NA,
  #                  'theo' = sum_theo,
  #                  'decomp.theo' = obj_teo,
  #                  'scales' = 2^seq(floor(log(length(data), 2))), 
  #                  'robust' = robust,
  #                  'eff' = eff,
  #                  'model.type' = model.type,
  #                  'compute.v' = compute.v,
  #                  'alpha' = alpha,
  #                  'expect.diff' = NA,
  #                  'N' = N,
  #                  'G' = G,
  #                  'H' = H,
  #                  'K' = K,
  #                  'model' = model,
  #                  'model.hat' = NA,
  #                  'starting' = model$starting,
  #                  'seed' = seed,
  #                  'freq' = freq,
  #                  'dr.slope' = NA), class = "gmwm")
  #   estimate = out[[1]]
  #   rownames(estimate) = model$process.desc
  #   colnames(estimate) = "Estimates" 
  #   
  # }else{
  
  # Standard GMWM using WV as input 
  wv_mat = cbind(wv$variance, wv$ci_low, wv$ci_high)
  out = .Call('_gmwm_gmwm_master_wv_cpp', PACKAGE = 'gmwm', wv_mat,
              wv$N, wv$mean_diff, Omega, wv$ranged,
              theta, desc, obj, model.type, starting = model$starting,
              p = alpha, compute_v = compute.v, K = K, H = H, G = G,
              robust = wv$robust, eff = wv$eff)
  
  estimate = out[[1]]
  rownames(estimate) = model$process.desc
  colnames(estimate) = "Estimates" 
  
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  colnames(init.guess) = "Starting" 
  #}
  
  # Convert from AR1 to GM
  if(detected_gm){
    estimate[,1] = conv.ar1.to.gm(estimate[,1], model$process.desc, freq)
    init.guess[,1] = conv.ar1.to.gm(init.guess[,1], model$process.desc, freq)
  }
  
  # Wrap this into the C++ Lib
  scales = .Call('_gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)
  
  # Create a new model object.
  model.hat = model
  
  model.hat$starting = F  
  
  model.hat$theta = as.numeric(estimate)
  
  # Release model
  #if(!all(model$process.desc %in% c('QN', 'WN', 'RW', 'DR'))){
  out = structure(list(estimate = estimate,
                       init.guess = init.guess,
                       wv.empir = out[[3]], 
                       ci.low = out[[4]], 
                       ci.high = out[[5]],
                       orgV = out[[7]],
                       V = out[[6]],
                       omega = out[[12]],
                       obj.fun = out[[11]],
                       theo = out[[9]],
                       decomp.theo = out[[10]],
                       scales = scales, 
                       robust = wv$robust,
                       eff = wv$eff,
                       model.type = model.type,
                       compute.v = compute.v,
                       alpha = alpha,
                       expect.diff = out[[8]],
                       N = N,
                       G = G,
                       H = H,
                       K = K,
                       model = model,
                       model.hat = model.hat,
                       starting = model$starting,
                       seed = seed,
                       freq = freq,
                       dr.slope = out[[13]],
                       wv = wv,
                       remove_scales = remove_scales), class = "gmwm")
  #}
  invisible(out)
}

#' Generalized Method of Wavelet Moments (GMWM) for IMUs, ARMA, SSM, and Robust
#' 
#' Performs estimation of time series models by using the GMWM estimator.
#' @param model      A \code{ts.model} object containing one of the allowed models.
#' @param data       A \code{matrix} or \code{data.frame} object with only column 
#'                   (e.g. \eqn{N \times 1}{ N x 1 }), a \code{lts} object,
#'                   or a \code{gts} object. 
#' @param model.type A \code{string} containing the type of GMWM needed:
#'                   \code{"imu"} or \code{"ssm"}.
#' @param compute.v  A \code{string} indicating the type of covariance matrix 
#'                   solver. Valid values are:
#'                    \code{"fast"}, \code{"bootstrap"}, 
#'                    \code{"diag"} (asymptotic diag), 
#'                    \code{"full"} (asymptotic full). By default, the program
#'                   will fit a "fast" model.
#' @param alpha      A \code{double} between 0 and 1 that correspondings to the
#'                   \eqn{\frac{\alpha}{2}}{alpha/2} value for the wavelet 
#'                   confidence intervals.
#' @param robust     A \code{boolean} indicating whether to use the robust
#'                   computation (\code{TRUE}) or not (\code{FALSE}).
#' @param eff        A \code{double} between 0 and 1 that indicates the
#'                   efficiency.
#' @param G          An \code{integer} to sample the space for IMU and SSM
#'                   models to ensure optimal identitability.
#' @param K          An \code{integer} that controls how many times the
#'                   bootstrapping procedure will be initiated.
#' @param H          An \code{integer} that indicates how many different
#'                   samples the bootstrap will be collect.
#' @param seed       An \code{integer} that controls the reproducibility of the
#'                   auto model selection phase.
#' @param freq       A \code{double} that indicates the sampling frequency. By
#'                   default, this is set to 1 and only is important if \code{GM()}
#'                   is in the model
#' @return A \code{gmwm} object with the structure: 
#' \describe{
#'  \item{estimate}{Estimated Parameters Values from the GMWM Procedure}
#'  \item{init.guess}{Initial Starting Values given to the Optimization Algorithm}
#'  \item{wv.empir}{The data's empirical wavelet variance}
#'  \item{ci.low}{Lower Confidence Interval}
#'  \item{ci.high}{Upper Confidence Interval}
#'  \item{orgV}{Original V matrix}
#'  \item{V}{Updated V matrix (if bootstrapped)}
#'  \item{omega}{The V matrix inversed}
#'  \item{obj.fun}{Value of the objective function at Estimated Parameter Values}
#'  \item{theo}{Summed Theoretical Wavelet Variance}
#'  \item{decomp.theo}{Decomposed Theoretical Wavelet Variance by Process}
#'  \item{scales}{Scales of the GMWM Object}
#'  \item{robust}{Indicates if parameter estimation was done under robust or classical}
#'  \item{eff}{Level of efficiency of robust estimation}
#'  \item{model.type}{Models being guessed}
#'  \item{compute.v}{Type of V matrix computation}
#'  \item{augmented}{Indicates moments have been augmented}
#'  \item{alpha}{Alpha level used to generate confidence intervals}
#'  \item{expect.diff}{Mean of the First Difference of the Signal}
#'  \item{N}{Length of the Signal}
#'  \item{G}{Number of Guesses Performed}
#'  \item{H}{Number of Bootstrap replications}
#'  \item{K}{Number of V matrix bootstraps}
#'  \item{model}{\code{ts.model} supplied to gmwm}
#'  \item{model.hat}{A new value of \code{ts.model} object supplied to gmwm}
#'  \item{starting}{Indicates whether the procedure used the initial guessing approach}
#'  \item{seed}{Randomization seed used to generate the guessing values}
#'  \item{freq}{Frequency of data}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' 
#' The V matrix is calculated by:
#' \eqn{diag\left[ {{{\left( {Hi - Lo} \right)}^2}} \right]}{diag[(Hi-Lo)^2]}.
#' 
#' The function is implemented in the following manner:
#' 1. Calculate MODWT of data with levels = floor(log2(data))
#' 2. Apply the brick.wall of the MODWT (e.g. remove boundary values)
#' 3. Compute the empirical wavelet variance (WV Empirical).
#' 4. Obtain the V matrix by squaring the difference of the WV Empirical's Chi-squared confidence interval (hi - lo)^2
#' 5. Optimize the values to obtain \eqn{\hat{\theta}}{theta^hat}
#' 6. If FAST = TRUE, return these results. Else, continue.
#' 
#'Loop  k = 1 to K
#' Loop h = 1 to H
#' 7. Simulate xt under \eqn{F_{\hat{\theta}}}{F_theta^hat}
#' 8. Compute WV Empirical
#' END
#' 9. Calculate the covariance matrix
#' 10. Optimize the values to obtain \eqn{\hat{\theta}}{theta^hat}
#'END
#' 11. Return optimized values.
#' 
#' 
#' The function estimates a variety of time series models. If type = "imu" or "ssm", then
#' parameter vector should indicate the characters of the models that compose the latent or state-space model. The model
#' options are:
#' \describe{
#'   \item{"AR1"}{a first order autoregressive process with parameters \eqn{(\phi,\sigma^2)}{phi, sigma^2}}
#'   \item{"GM"}{a guass-markov process \eqn{(\beta,\sigma_{gm}^2)}{beta, sigma[gm]^2}}
#'   \item{"ARMA"}{an autoregressive moving average process with parameters \eqn{(\phi _p, \theta _q, \sigma^2)}{phi[p], theta[q], sigma^2}}
#'   \item{"DR"}{a drift with parameter \eqn{\omega}{omega}}
#'   \item{"QN"}{a quantization noise process with parameter \eqn{Q}}
#'   \item{"RW"}{a random walk process with parameter \eqn{\sigma^2}{sigma^2}}
#'   \item{"WN"}{a white noise process with parameter \eqn{\sigma^2}{sigma^2}}
#' }
#' If only an ARMA() term is supplied, then the function takes conditional least squares as starting values
#' If robust = TRUE the function takes the robust estimate of the wavelet variance to be used in the GMWM estimation procedure.
#' @export
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' data = gen_gts(n, AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 1))
#' 
#' # Models can contain specific parameters e.g.
#' adv.model = gmwm(AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 0.01),
#'                             data)
#' 
#' # Or we can guess the parameters:
#' guided.model = gmwm(AR1() + WN(), data) 
#' 
#' # Want to try different models? 
#' guided.ar1 = gmwm(AR1(), data)
#' 
#' # Faster:
#' guided.ar1.wn.prev = update(guided.ar1, AR1()+WN())
#' 
#' # OR 
#' 
#' # Create new GMWM object. 
#' # Note this is SLOWER since the Covariance Matrix is recalculated.
#' guided.ar1.wn.new = gmwm(AR1()+WN(), data)
#'  
#' # ARMA case
#' set.seed(1336)
#' data = gen_gts(n, ARMA(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488),
#'               sigma2 = 0.1796))
#' #guided.arma = gmwm(ARMA(2,2), data, model.type="ssm")
#' adv.arma = gmwm(ARMA(ar=c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), sigma2=0.1796),
#'                 data, model.type="ssm")
gmwm_data = function(model, data, model.type = "imu", compute.v = "auto", remove_scales = NULL, 
                     Omega = NULL, robust = FALSE, eff = 0.6, alpha = 0.05, seed = 1337, 
                     G = NULL, K = 1, H = 100, freq = 1){
  
    wv = wvar(data, robust = robust, eff = eff)
    
    gmwm_wvar(model = model, 
              wv = wv, 
              model.type = model.type, 
              compute.v = compute.v, 
              remove_scales = remove_scales, 
              Omega = Omega,
              alpha = alpha, 
              seed = seed, 
              G = G, 
              K = K, 
              H = H,
              freq = freq)
}


# TO DOCUMENT + USEMETHOD (SEE WVAR)
#' @export
mgmwm = function(model, input, model.type = "imu", compute.v = "auto", remove_scales = NULL, trim = FALSE,
                 Omega = NULL, alpha = 0.05, seed = 1337, 
                 G = NULL, K = 1, H = 100, freq = 1){
  if (class(input[[1]])[1] == "wvar"){
    # Number of samples
    M = length(input)
    
    # Find largest J
    Js =  rep(NA, M)
    for (i in 1:M){
      Js[i] = input[[i]]$J
    }
    Jmax = max(Js)
    
    # Construct WV array
    wv_array = array(NA, c(Jmax, 3, M))
    
    for (i in 1:M){
      wv_array[1:Js[i], 1:3, i] = cbind(input[[i]]$variance, input[[i]]$ci_low, input[[i]]$ci_high)
    }
    
    # Compute Omega
    if (is.null(Omega)){
      Omega = matrix(0, Jmax, Jmax)
      for (i in 1:M){
        Omega[1:Js[i], 1:Js[i]] = Omega[1:Js[i], 1:Js[i]] + diag(1/(input[[i]]$ci_high - input[[i]]$ci_low)^2)
      }
    }
    
    # Compute mean WV
    wv_input = input[[which.max(Js)]]
    wv_input$variance = apply(wv_array[,1,], 1, mean, na.rm = TRUE, trim = as.integer(trim))
    
    # Fit Standard GMWM
    fit = gmwm_wvar(model = model, 
              wv = wv_input, 
              model.type = model.type, 
              compute.v = compute.v, 
              remove_scales = remove_scales, 
              Omega = Omega,
              alpha = alpha, 
              seed = seed, 
              G = G, 
              K = K, 
              H = H,
              freq = freq)
    
    out = list(wv = wv_input, wv_array = wv_array, fit = fit, names = names(input))
    class(out) = "mgmwm"
    return(out)
  }else{
    print("Not implemented :(")
  }
}


# TO DOCUMENT + USEMETHOD (SEE WVAR)
#' @export
gmwm = function(model, input, model.type = "imu", compute.v = "auto", remove_scales = NULL, 
                     Omega = NULL, robust = FALSE, eff = 0.6, alpha = 0.05, seed = 1337, 
                     G = NULL, K = 1, H = 100, freq = 1){
  
  # ADD CHECKS
  if (class(input)[1] == "wvar"){
    return(gmwm_wvar(model = model, 
              wv = input, 
              model.type = model.type, 
              compute.v = compute.v, 
              remove_scales = remove_scales, 
              Omega = Omega,
              alpha = alpha, 
              seed = seed, 
              G = G, 
              K = K, 
              H = H,
              freq = freq))
  }else{
    return(gmwm_data(model = model, data = input, model.type = model.type, 
              compute.v = compute.v, remove_scales = remove_scales, 
              Omega = Omega, robust = robust, eff = eff, alpha = alpha, 
              seed = seed, G = G, K = K, H = H, freq = freq))
  }
}


#' Print gmwm object
#'
#' Displays information about GMWM object
#' @method print gmwm
#' @export
#' @keywords internal
#' @param x   A \code{GMWM} object
#' @param ... Other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' \dontrun{
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen_gts(n, AR1(phi=.1, sigma2 = 1) + AR1(phi=0.95, sigma2 = .1))
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' print(mod)
#' }
print.gmwm = function(x, ...){
  cat("Model Information: \n")
  print(x$estimate)

  cat("\n* The initial values of the parameters used in the minimization of the GMWM objective function \n  were", 
      {if(x$starting) paste0("generated by the program underneath seed: ",x$seed,".") else "supplied by you."},"\n\n")
}


#' Summary of GMWM object
#'
#' Displays summary information about GMWM object
#' @method summary gmwm
#' @export
#' @param object       A \code{GMWM} object
#' @param inference    A value containing either: NULL (auto), TRUE, or FALSE
#' @param bs.gof       A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.gof.p.ci  A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.theta.est A value containing either: NULL (auto), TRUE, FALSE
#' @param bs.ci        A value containing either: NULL (auto), TRUE, FALSE
#' @param B            An \code{int} that indicates how many bootstraps should be performed.
#' @param ...          Other arguments passed to specific methods
#' @return A \code{summary.gmwm} object with:
#' \describe{
#'  \item{estimate}{Estimated Theta Values}
#'  \item{testinfo}{Goodness of Fit Information}
#'  \item{inference}{Inference performed? T/F}
#'  \item{bs.gof}{Bootstrap GOF? T/F}
#'  \item{bs.gof.p.ci}{Bootstrap GOF P-Value CI? T/F}
#'  \item{bs.theta.est}{Bootstrap Theta Estimates? T/F}
#'  \item{bs.ci}{Bootstrap CI? T/F}
#'  \item{starting}{Indicates if program supplied initial starting values}
#'  \item{seed}{Seed used during guessing / bootstrapping}
#'  \item{obj.fun}{Value of obj.fun at minimized theta}
#'  \item{N}{Length of Time Series}
#' }
#' @author JJB
#' @examples
#' \dontrun{
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen_gts(n, AR1(phi=.1, sigma2 = 1) + AR1(phi=0.95, sigma2 = .1))
#' mod = gmwm(AR1()+AR1(), data = xt, model.type = "imu")
#' summary(mod)
#' }
summary.gmwm = function(object, inference = NULL,  
                        bs.gof = NULL,  bs.gof.p.ci = NULL, 
                        bs.theta.est = NULL, bs.ci = NULL,
                        B = 100, ...){
  
  # Set a different seed to avoid dependency.
  set.seed(object$seed+5)
  
  out = object$estimate
  
  N = object$N
  
  # Enable values if small time series.
  auto = if(N > 10000) FALSE else TRUE
  
  # Auto set values
  if(is.null(inference)){
    inference = auto
  }

  if(is.null(bs.gof)){
    bs.gof= if(inference) auto else F
  }
  
  if(is.null(bs.gof.p.ci)){
    bs.gof.p.ci = if(inference) auto else F
  }
  
  if(is.null(bs.theta.est)){
    bs.theta.est = if(inference) auto else F
  }
  
  if(is.null(bs.ci)){
    bs.ci = if(inference) auto else F
  }
  
  if("ARMA" %in% object$model$desc){
    if(bs.ci == FALSE){
      warning(paste0("The numerical derivative of ARMA(p,q), where p > 1 and q > 1, may be inaccurate leading to inappropriate CIs.\n",
              "Consider using the bs.ci = T option on the summary function."))
    }
  }
  
  if(inference){
    
    # Convert from GM to AR1
    if(any(object$model$desc == "GM")){
      object$estimate[,1] = conv.gm.to.ar1(object$estimate[,1], object$model$process.desc, object$freq)
    }
    
    mm = .Call('_gmwm_get_summary', PACKAGE = 'gmwm',object$estimate,
                                                    object$model$desc, object$model$obj.desc,
                                                    object$model.type, 
                                                    object$wv.empir, object$theo,object$scales,
                                                    object$V, solve(object$orgV), object$obj.fun,
                                                    N, object$alpha,
                                                    object$robust, object$eff,
                                                    inference, F, # fullV is always false. Need same logic updates.
                                                    bs.gof, bs.gof.p.ci, bs.theta.est, bs.ci,
                                                    B)
  }else{
    mm = vector('list',3)
    mm[1:3] = NA
  }
  
  if(inference){
    out.coln = colnames(out)
    out = cbind(out, mm[[1]])
    colnames(out) = c(out.coln, "CI Low", "CI High", "SE")
    # Convert from AR1 to GM
    idx_gm = (object$model$desc == "GM")
    if(any(idx_gm)) {
      out[,2:3] = apply(out[,2:3], 2, FUN = conv.ar1.to.gm,
                                      process.desc = object$model$process.desc, 
                                      freq = object$freq)
      # To do: Add delta method transform here for sigma2
    }
    
  }
  
  x = structure(list(estimate=out, 
                     testinfo=mm[[2]],
                     inference = inference, 
                     bs.gof = bs.gof,
                     bs.gof.p.ci = bs.gof.p.ci,
                     bs.theta.est = bs.theta.est, 
                     bs.ci = bs.ci,
                     starting = object$starting,
                     seed = object$seed,
                     obj.fun = object$obj.fun,
                     N = N,
                     freq = object$freq), class = "summary.gmwm")
    
  x
}

#' Print summary.gmwm object
#'
#' Displays summary information about GMWM object
#' @method print summary.gmwm
#' @export
#' @keywords internal
#' @param x   A \code{GMWM} object
#' @param ... Other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' \dontrun{
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen_gts(n, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' mod = gmwm(AR1() + AR1(), data = xt, model.type = "imu")
#' summary(mod)
#' }
print.summary.gmwm = function(x, ...){
  
  cat("Model Information: \n")
  print(x$estimate)
  if(x$bs.theta.est){
    cat("\n> The parameter estimates shown are bootstrapped! To use these results, please save the summary object.")
  }
  
  cat("\n* The initial values of the parameters used in the minimization of the GMWM objective function \n  were", 
      {if(x$starting) paste0("generated by the program underneath seed: ",x$seed,".") else "given by YOU!"},"\n\n")

  cat(paste0("Objective Function: ", round(x$obj.fun,4),"\n\n"))
  
    
  if(x$inference){
    cat(paste0({if(x$bs.gof) "Bootstrapped" else "Asymptotic"}," Goodness of Fit: \n"))
    if(x$bs.gof){
      cat(paste0("Test Statistic: ", round(x$obj.fun,2),"\n",
                 "P-Value: ", round(x$testinfo[1],4)), 
                {if(x$bs.gof.p.ci) paste0(" CI: (", round(x$testinfo[2],4),", ", round(x$testinfo[3],4), ")")})

    }else{
      cat(paste0("Test Statistic: ", round(x$testinfo[1],2),
          " on ",x$testinfo[3]," degrees of freedom\n",
          "The resulting p-value is: ", round(x$testinfo[2],4)))
    }
    cat("\n\n")
  }
  
  if(x$bs.gof || x$bs.theta.est)
   cat(paste0("\nTo replicate the results, use seed: ",x$seed, "\n"))
}

#' Predict future points in the time series using the solution of the
#' Generalized Method of Wavelet Moments
#' 
#' Creates a prediction using the estimated values of GMWM through the 
#' ARIMA function within R.
#' @param object       A \code{\link{gmwm}} object 
#' @param data.in.gmwm The data SAME EXACT DATA used in the GMWM estimation
#' @param n.ahead      Number of observations to forecast
#' @param ...          Additional parameters passed to ARIMA Predict
#' @return A \code{predict.gmwm} object with:
#' \describe{
#' \item{pred}{Predictions}
#' \item{se}{Standard Errors}
#' \item{resid}{Residuals from ARIMA ML Fit}
#' }
#' @seealso \code{\link{gmwm}}, \code{\link{ARMA}}
#' @export
#' @examples 
#' # Simulate an ARMA Process
#' xt = gen_gts(1000, ARMA(ar=0.3, ma=0.6, sigma2=1))
#' model = gmwm(ARMA(1,1), xt)
#' 
#' # Make prediction
#' predict(model, xt, n.ahead = 1)
predict.gmwm = function(object, data.in.gmwm, n.ahead = 1, ...){
  
  ts.mod = object$model
  
  if(length(ts.mod$desc) > 1 || ts.mod$desc != "SARIMA")
    stop("The predict function only works with stand-alone SARIMA models.")
  
  objdesc = ts.mod$obj.desc[[1]]

  # Unpack ts object
  p = objdesc[1]
  q = objdesc[2]
  P = objdesc[3]
  Q = objdesc[4]
  s = objdesc[6] # Set to 0 (handled in ARIMA)
  d = objdesc[7]
  D = objdesc[8]
  
  # Make an ARIMA object
  mod = arima(data.in.gmwm, order = c(p, d, q),
              list(order = c(P, D, Q), period = s),
              method = "ML",
              fixed = object$estimate[1:(p+q+P+Q)],
              transform.pars = F,
              include.mean = F)
  
  # Predict off of ARIMA
  pred = predict(mod, n.ahead = n.ahead, newxreg = NULL,
                 se.fit = TRUE, ...)
  
  # Format Results
  structure(list(pred = pred$pred,
                   se = pred$se,
                resid = mod$residuals),
            class = "predict.gmwm")
                          
}

# TO DOCUMENT + ADD OPTIONS
#' @title plot method for GMWM
#' @method  plot gmwm
#' @export
plot.gmwm = function(x, decomp = TRUE,
                     add_legend_gwmw = TRUE, 
                     ylab_gmwm = NULL,
                     col_theo = "#F47F24",
                     pch_theo = 1,
                     cex_theo = 1.25,
                     lwd_theo = 1.25,
                     add_legend = TRUE,
                     legend_pos = NULL,
                     col_wv = NULL,
                     col_ci = NULL){
  
  
  # SOME CHECKS SHOULD BE ADDED
  
  if (is.null(ylab_gmwm)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }else{
    ylab = ylab_gmwm
  }
  
  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }
  
  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }
  
  par(mar=c(5.1, 5.1, 1, 12))
  plot(x$wv, legend_position = NA, ylab = ylab, 
       col_ci = col_ci, col_wv = col_wv)
  
  if (decomp){
    M = ncol(x$decomp.theo)
    wv_cols = ggColor(M + 1)[-1]
    
    for (m in 1:M){
      lines(x$scales, x$decomp.theo[,m], 
            col = wv_cols[m])
    }
  }
  
  lines(x$scales, apply(x$decomp.theo, 1, sum), 
        col = col_theo, pch = pch_theo, type = "b", 
        cex = cex_theo, lwd = lwd_theo)
  
  if (!is.null(x$remove_scales)){
    # ADD CHECKS HERE
    points(x$scales[x$remove_scales], x$wv.empir[x$remove_scales], col = "red3", pch = 4, lwd = 3, cex = 1.5)
  }
  
  # TO DO: CHECK PLACEMEMENT OF LEGEND
  # Legend 
  if (add_legend == TRUE){
    
    if (x$robust == TRUE){
      legend_names = c("Robust GMWM:", "Implied WV")
    }else{
      legend_names = c("Standard GMWM:", "Implied WV")
    }
    
    
    col_leg = c(NA, col_theo)
    pch_leg = c(NA, pch_theo)
    cex_leg = c(NA, 1.2*cex_theo)
    lwd_leg = c(NA, lwd_theo)
    lty_leg = c(NA, 1)
    pt_lwd = c(1, 1)
    
    CI_conf = 1 - x$alpha
    
    if (x$robust == TRUE){
      wv_title_part1 = "Empirical Robust WV "
    }else{
      wv_title_part1 = "Empirical WV "
    }
    
    legend_names = c(legend_names, 
                     as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))), 
                     as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")"))))
    
    col_leg = c(col_leg, col_wv, col_ci) 
    lty_leg = c(lty_leg, 1, NA)
    pch_leg = c(pch_leg, 16, 15)
    cex_leg = c(cex_leg, 1.25, 3)
    lwd_leg = c(lwd_leg, 1, NA) 
    pt_lwd = c(pt_lwd, 1, 1)
    
    if (decomp){
      legend_names = c(legend_names, NA, "Latent Processes:", strsplit(x$model$print,split=' + ', fixed=TRUE)[[1]])
      col_leg = c(col_leg, NA, NA, wv_cols)
      lty_leg = c(lty_leg, NA, NA, rep(1, M))
      pch_leg = c(pch_leg, NA, NA, rep(NA, M))
      cex_leg = c(cex_leg, NA, NA, rep(1, M))
      lwd_leg = c(lwd_leg, NA, NA, rep(1, M))
      pt_lwd = c(pt_lwd, NA, NA, rep(1, M))
    }
    
    if (!is.null(x$remove_scales)){
      legend_names = c(legend_names, NA, "Discarded Scales")
      col_leg = c(col_leg, NA, "red3")
      lty_leg = c(lty_leg, NA, NA)
      pch_leg = c(pch_leg, NA, 4)
      cex_leg = c(cex_leg, NA, 1.5)
      lwd_leg = c(lwd_leg, NA, NA)
      pt_lwd = c(pt_lwd, NA, 3)
    }
    
    legend("topleft", legend_names,  inset=c(1.03,0.1),
           bty='n', xpd=TRUE,
           col = col_leg,
           pch = pch_leg,
           pt.cex = cex_leg,
           pt.lwd = pt_lwd,
           lwd = lwd_leg,
           cex = 0.95)
  }
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}


# TO DOCUMENT + ADD OPTIONS
#' @title plot method for MGMWM
#' @method  plot mgmwm
#' @export
plot.mgmwm = function(x, decomp = TRUE,
                     add_legend_gwmw = TRUE, 
                     ylab_gmwm = NULL,
                     col_theo = "#F47F24",
                     pch_theo = 1,
                     cex_theo = 1.25,
                     lwd_theo = 1.25,
                     add_legend = TRUE,
                     legend_pos = NULL,
                     col_wv = NULL,
                     col_ci = NULL){
  
  
  # SOME CHECKS SHOULD BE ADDED
  
  if (is.null(ylab_gmwm)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }else{
    ylab = ylab_gmwm
  }
  
  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }
  
  par(mar=c(5.1, 5.1, 1, 12))
  plot(x$wv, legend_position = NA, ylab = ylab, 
       col_wv = col_wv, ci_wv = FALSE)
  
  if (decomp){
    M1 = ncol(x$fit$decomp.theo)
    M2 = dim(x$wv_array)[3]
  }else{
    M1 = 0
    M2 = dim(x$wv_array)[3]
  }
  
  M = M1 + M2
  wv_cols = ggColor(M)
  wv_cols_trans = ggColor(M, alpha = 0.15)
  
  scales = 2^(1:dim(x$wv_array)[1])/x$fit$freq
  
  for (i in 1:M2){
    wv_current = na.omit(x$wv_array[,,i])
    J_current = nrow(wv_current)
    lines(scales[1:J_current], wv_current[,1], col = wv_cols[i], lwd = 1.25, lty = 2)
    polygon(c(scales[1:J_current], rev(scales[1:J_current])),
            c(wv_current[,2], rev(wv_current[,3])), col = wv_cols_trans[i],
            border = NA)
  }
    
  for (i in 1:M1){
      lines(scales, x$fit$decomp.theo[,i], 
            col = wv_cols[i + M2], lwd = 1.5)
  }
  
  lines(scales, apply(x$fit$decomp.theo, 1, sum), 
        col = col_theo, pch = pch_theo, type = "b", 
        cex = cex_theo, lwd = lwd_theo)
  
  if (!is.null(x$fit$remove_scales)){
    # ADD CHECKS HERE
    points(x$fit$scales[x$fit$remove_scales], x$fit$wv.empir[x$fit$remove_scales], col = "red3", pch = 4, lwd = 3, cex = 1.5)
  }
  
  # TO DO: CHECK PLACEMEMENT OF LEGEND
  # Legend 
  if (add_legend == TRUE){
    
    if (x$fit$robust == TRUE){
      legend_names = c("Robust MGMWM:", "Implied WV")
    }else{
      legend_names = c("Standard MGMWM:", "Implied WV")
    }
    
    
    col_leg = c(NA, col_theo)
    pch_leg = c(NA, pch_theo)
    cex_leg = c(NA, 1.2*cex_theo)
    lwd_leg = c(NA, lwd_theo)
    lty_leg = c(NA, 1)
    pt_lwd = c(1, 1)
    
    CI_conf = 1 - x$alpha
    
    if (x$fit$robust == TRUE){
      wv_title_part1 = "Mean Robust WV "
    }else{
      wv_title_part1 = "Mean WV "
    }
    
    legend_names = c(legend_names, 
                     as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))))
    
    col_leg = c(col_leg, col_wv) 
    lty_leg = c(lty_leg, 1)
    pch_leg = c(pch_leg, 16)
    cex_leg = c(cex_leg, 1.25)
    lwd_leg = c(lwd_leg, 1) 
    pt_lwd = c(pt_lwd, 1)
    
    if (decomp){
      legend_names = c(legend_names, NA, "Latent Processes:", strsplit(x$fit$model$print,split=' + ', fixed=TRUE)[[1]])
      col_leg = c(col_leg, NA, NA, wv_cols[(M2 + 1):(M1 + M2)])
      lty_leg = c(lty_leg, NA, NA, rep(1, M1))
      pch_leg = c(pch_leg, NA, NA, rep(NA, M1))
      cex_leg = c(cex_leg, NA, NA, rep(1, M1))
      lwd_leg = c(lwd_leg, NA, NA, rep(1, M1))
      pt_lwd = c(pt_lwd, NA, NA, rep(1, M1))
    }
    
    if (!is.null(x$fit$remove_scales)){
      legend_names = c(legend_names, NA, "Discarded Scales")
      col_leg = c(col_leg, NA, "red3")
      lty_leg = c(lty_leg, NA, NA)
      pch_leg = c(pch_leg, NA, 4)
      cex_leg = c(cex_leg, NA, 1.5)
      lwd_leg = c(lwd_leg, NA, NA)
      pt_lwd = c(pt_lwd, NA, 3)
    }
    
    if (is.null(x$names)){
      x$names = as.character(1:M2)
    }
    legend_names = c(legend_names, NA, "Samples:", x$names)
    col_leg = c(col_leg, NA, NA, wv_cols_trans[(1):(M2)])
    lty_leg = c(lty_leg, NA, NA, rep(2, M2))
    pch_leg = c(pch_leg, NA, NA, rep(15, M2))
    cex_leg = c(cex_leg, NA, NA, rep(2, M2))
    lwd_leg = c(lwd_leg, NA, NA, rep(1, M2))
    pt_lwd = c(pt_lwd, NA, NA, rep(1, M2))
    
    legend("topleft", legend_names,  inset=c(1.03,0.1),
           bty='n', xpd=TRUE,
           col = col_leg,
           pch = pch_leg,
           pt.cex = cex_leg,
           pt.lwd = pt_lwd,
           lwd = lwd_leg,
           cex = 0.95)
  }
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}



#' Graphically Compare GMWM Model Fit
#'
#' Creates a table of graphs to compare GMWM model fit.
#' @return A ggplot2 panel containing the graphs of gmwm objects.
#' @param ... Several \code{gmwm} objects.
#' @param display.model A \code{boolean} indicating whether the model should be displayed in the facet label.
#' @param show.theo.wv A \code{boolean} indicating whether the theoretical WV should be plotted in the lower triangular area. See \code{details}.
#' @param facet.label  A \code{character vector} indicating what should be displayed in the facet labels.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of confidence interval.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of the confidence interval.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @export
#' @author Wenchao
#' @details 
#' This function has been updated to deal with any \code{gmwm} object. The old contraint that supplied
#' objects must be constrcuted by the same data is no longer needed.
#' 
#' The parameter \code{show.theo.wv} will be automatically set to TRUE, if the supplied \code{gmwm} objects
#' are constructed by the same data. This aims to mimic the behaviour in the old version of \code{compare_models}.
#' 
#' The default aesthetical setting is designed to work well with 2 and 3 objects. Users are expected to change
#' the color/point/line settings if they want to supply more than 3 objects.
#' 
#' If two models have the same attribute values, for example, same empirical WV, then their aethetics will
#' be set to the same, i.e. same line type, line color, point size, point shape, etc. 
#' 
#' @examples
#' \dontrun{
#' #AR
#' set.seed(8836)
#' x1 = gen_gts(1000, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' x2 = gen_gts(2000, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' 
#' GMWM1 = gmwm(AR1(), data = x1)
#' GMWM2 = gmwm(2*AR1(), data = x2)
#' 
#' compare_models(GMWM1, GMWM2, show.theo.wv = T, transparence = 0.2, 
#'                facet.label = c('model1', 'model2'))
#'                
#' compare_models(GMWM1, GMWM2, CI.color = c('black','red'), 
#'                point.size = c(3, 3, 0, 0, 0, 0, 2, 2),
#'                line.color = c('black', 'red', 'grey', 'pink',
#'                               'grey', 'pink', 'purple', 'green')) 
#' #in the order of emp. WV, lower bound, higher bound, theo. WV
#' }
compare_models = function(..., display.model = T, show.theo.wv = F, facet.label = NULL, 
                        background = 'white', transparence = 0.05, CI.color = NULL,
                        line.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                        title = "Comparison of GMWM Models", title.size= 18, 
                        axis.label.size = 16, axis.tick.size = 11, 
                        facet.label.size = 13, facet.label.background = "#003C7D33",
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu))){
  scales=variable=l_value=h_value=low=.x=NULL
  
  # S1: Checking statement (Reset it to default setting if user passes wrong values)
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...) 
  numObj = length(obj_list)
  object.names = as.character(substitute(...()))
  
  if(numObj<1){
    stop('At least one model must be supplied.')
  }
  
  if(numObj == 1){
    #plot.gmwm will do the parameter checking
    #other parameters are not listed here, and they cannot be passed to plot.gmww by '...'
    warning('One object is supplied. You are actually calling plot.gmwm().')
    
    if(is.null(CI.color)) {CI.color="#003C7D"}
    
    plot.gmwm(x = obj_list[[1]], process.decomp = F, background = background, transparence = transparence, 
              CI.color = CI.color, CI = T, bw = F, line.type = line.type,
              line.color = line.color, point.size = point.size, point.shape = point.shape, title = title,
              title.size = title.size, axis.label.size = axis.label.size, axis.tick.size = axis.tick.size,
              axis.x.label = axis.x.label, axis.y.label = axis.y.label)
    
  }else{
    # check whether supplied objects are valid gmwm objects; if not, stop
    sapply(obj_list, FUN = function(x){ 
      if( !is.gmwm(x) ){
        stop("The function can only work on 'gmwm' object.")} })
    
    # freq conversion
    for (i in 1:numObj ){
      obj_list[[i]]$scales = obj_list[[i]]$scales/obj_list[[i]]$freq
    }
    
    # check parameter
    params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape', 'facet.label')
    requireLength = c(4*numObj, 4*numObj, numObj, 4*numObj, 4*numObj, numObj)
    default = list(NULL, NULL, NULL, NULL, NULL, NULL)
    nullIsFine = c(rep(T,6))
    
    checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
    
    # S2: Auto-select parameters, if not provided by users
    
    # decide what should appear in facet label
    if(is.null(facet.label) && display.model){
      # melt function cannot deal with expression
      object.names = sapply(obj_list, FUN = function(x){as.character( getModel.gmwm(x) ) } )
    }else if(!is.null(facet.label)){
      object.names = facet.label 
    }
    
    # Problem: When object names are the same, function will not work
    # Solution: Add space after object names
    object.names = addSpaceIfDuplicate(object.names)
    
    # in the order: WV low high theo
    if( is.null(line.type) ){line.type = c(rep('solid',numObj),#WV
                                           rep('dotted', numObj),#low
                                           rep('dotted', numObj),#high
                                           rep('solid', numObj))}#theo
    
    if(is.null(CI.color)){
      CI.color = switch(as.character(numObj),
                        '2' = c("#003C7D","#F47F24"), #UIUC
                        '3' = c("#003C7D","#F47F24","#34A853"), 
                        ggColor(numObj) )
      
    }
    
    Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
    if(is.null(line.color)) {
      if (numObj <= 3) {
        line.color = switch(as.character(numObj),
                            '2' = c(rep(CI.color,3), c("#E41A1C","#4DAF4A")),#set 1
                            '3' = c(rep(CI.color,3), c("#E41A1C", "#4DAF4A", "#341456"))) #set 1 and purple
      }else{
        modulus = numObj %/% 8
        remainder = numObj %% 8
        theo.color =  c(rep(Set1, times = modulus), Set1[1:remainder])
        line.color = c(rep(CI.color, 3), #WV,low,high
                       theo.color) #theo
        if(numObj>8){
          warning('More than 8 objects are supplied, but the palette has only 8 colors')
        }
        
      }
    }
    
    #point.size will auto-decrease as the number of models increases
    if(is.null(point.size)){
      
      if(numObj > 5){temp.point.size = 1
      }else{
        temp.point.size = 3-0.5*(numObj-1)}
      
      point.size = c(rep(temp.point.size, numObj),#WV
                     rep(0, numObj),#low
                     rep(0, numObj),#high
                     rep(temp.point.size, numObj))#theo
      
    }
    
    if(is.null(point.shape)){point.shape = c(rep(20, numObj),#WV
                                             rep(46, numObj),#low
                                             rep(46, numObj),#high
                                             rep(1, numObj))}#theo
    
    # Check whether gmwm object have same attributes. If yes, then set the aethetics (color, shape etc.) the same.
    equal.attr = checkEquality(obj_list)
    # turn show.theo.wv on if gmwm objects are constructed by same data, i.e. same scales, same wv.empir, same bounds
    same.data = rbind(equal.attr$scales, equal.attr$wv.empir, equal.attr$bounds)
    if(all(same.data == T)){
      if(!show.theo.wv){warning("'show.theo.wv' is automatically set to TRUE.")}
      show.theo.wv = T
    }
    
    target.attrs = c('wv.empir', 'bounds', 'bounds', 'theo')
    
    for(ini.aes in c('line.type', 'line.color', 'point.size', 'point.shape')){
      for(counter in 0:(4*numObj-1)){
        final.aes = get(ini.aes)
        modulus2 = counter%/%numObj
        remainder2 = counter%%numObj + 1
        
        target.attr = switch(as.character(modulus2),
                             '0' = 'wv.empir',
                             '1' = 'bounds',
                             '2' = 'bounds',
                             '3' = 'theo')
        
        if(remainder2>1){
          
          for(i in (remainder2-1):1){ #only compare with previous gmwm objects
            if(equal.attr[[target.attr]][remainder2, i]){
              final.aes[counter+1] = final.aes[numObj*modulus2 + i]
            }
          }
        }
        
        assign(ini.aes, final.aes)
        
      }#end of 'counter' loop
    }#end of 'ini.aes' loop
    
    # CI.color still needs to be checked
    for(i in 2:numObj){
      for(j in 1:(i-1)){#only compare with previous gmwm objects
        
        if(equal.attr$bounds[i, j]){
          CI.color[i] = CI.color[j]
          
        }
      }
    }
    
    # S3: Rearrange the data into a data frame which can be passed to next step
    each.len = sapply(X = obj_list, FUN = function(x) length(x$wv.empir) )
    max.len = max(each.len)
    total.len = max.len*numObj^2
    
    #Choose the column names for empty data frame
    nlen = nchar(numObj)
    scales_names = sprintf(paste0("scales%0",nlen,"d"),1:numObj)
    WV_names = sprintf(paste0("WV%0",nlen,"d"),1:numObj)
    low_names = sprintf(paste0("low%0",nlen,"d"),1:numObj)
    high_names = sprintf(paste0("high%0",nlen,"d"),1:numObj)
    z_theo_names = sprintf(paste0("z_theo%0",nlen,"d"),1:numObj)
    col.names = c(scales_names, 'dataset_h', 'dataset_v', 
                  WV_names, low_names, high_names, z_theo_names)
    
    # Initialize empty data frame with right number of rows
    num.col = 2 + 5*numObj
    obj = as.data.frame(matrix(NA, ncol = num.col, nrow = total.len), stringsAsFactors = FALSE)
    colnames(obj) = col.names
    
    for(i in 1:numObj){
      obj[,scales_names[i]] = rep( autofill(obj_list[[i]]$scales, max.len), times = numObj^2 )
    }
    
    #put data into data frame
    t = 1
    d = max.len
    
    for (i in 1:numObj){
      for(j in 1:numObj){
        
        ## Diagonal ========================
        if(i == j){
          #compare to itself
          #set some values to NA
          
          WV_name = WV_names[i]
          low_name = low_names[i] 
          high_name = high_names[i] 
          z_theo_name = z_theo_names[i]
          
          obj[t:(t+d-1),c('dataset_h', 'dataset_v', 
                          WV_name, low_name, high_name, z_theo_name)] = 
            
            data.frame(object.names[i],
                       object.names[j],
                       
                       autofill(obj_list[[i]]$wv.empir, max.len),
                       autofill(obj_list[[i]]$ci.low, max.len),
                       autofill(obj_list[[i]]$ci.high, max.len),
                       autofill(obj_list[[i]]$theo, max.len),
                       
                       stringsAsFactors = F)
          t = t +d
          
        }else if (i > j ){
          ## lower triagular ========================
          
          WV_name = WV_names[c(i,j)]
          low_name = low_names[c(i,j)] 
          high_name = high_names[c(i,j)] 
          
          # WV
          WV.i = autofill( obj_list[[i]]$wv.empir, max.len)
          if(equal.attr$wv.empir[i,j]){
            WV.j = NA
          }else{
            WV.j = autofill( obj_list[[j]]$wv.empir, max.len)
          }
          
          # bound
          low.i = autofill( obj_list[[i]]$ci.low, max.len)
          high.i = autofill( obj_list[[i]]$ci.high, max.len)
          
          if(equal.attr$bounds[i,j]){
            low.j = NA
            high.j = NA
          }else{
            low.j = autofill( obj_list[[j]]$ci.low, max.len)
            high.j = autofill( obj_list[[j]]$ci.high, max.len)
          }
          
          # show.theo.wv = T
          if(show.theo.wv){
            z_theo_name = z_theo_names[c(i,j)]
            
            theo.i = autofill( obj_list[[i]]$theo, max.len )
            if(equal.attr$theo[i,j]){
              theo.j = NA
            }else{
              theo.j = autofill( obj_list[[j]]$theo, max.len)
            }
            
            obj[t:(t+d-1), c('dataset_h', 'dataset_v', WV_name, low_name, high_name, z_theo_name)] = 
              
              data.frame(object.names[i],
                         object.names[j],
                         
                         WV.i,
                         WV.j,
                         
                         low.i,
                         low.j,
                         
                         high.i,
                         high.j,
                         
                         theo.i,
                         theo.j, stringsAsFactors = F)
            
          }else{
            # show.theo.wv = F
            
            obj[t:(t+d-1),c('dataset_h', 'dataset_v', WV_name, low_name, high_name)] = 
              
              data.frame(object.names[i],
                         object.names[j],
                         
                         WV.i,
                         WV.j,
                         
                         low.i,
                         low.j,
                         
                         high.i,
                         high.j, stringsAsFactors = F)
          }
          
          t = t +d
        }else{
          ## upper triagular ======================== 
          z_theo_name = z_theo_names[c(i,j)]
          
          theo.i = autofill( obj_list[[i]]$theo, max.len )
          if(equal.attr$theo[i,j]){
            theo.j = NA
          }else{
            theo.j = autofill( obj_list[[j]]$theo, max.len )
          }
          
          obj[t:(t+d-1), c('dataset_h','dataset_v', z_theo_name)] = 
            
            data.frame(object.names[i],
                       object.names[j],
                       
                       theo.i,
                       theo.j, stringsAsFactors = F)
          
          t = t +d
        }
      } # end of j loop
    }# end of i loop
    
    # change the order of facet label
    obj$dataset_h = factor(obj$dataset_h, levels = object.names )
    obj$dataset_v = factor(obj$dataset_v, levels = object.names )
    
    # 1. data frame that is used to plot lines
    line_df.temp = melt(obj, id.vars = c(scales_names, 'dataset_v', 'dataset_h'), factorsAsStrings = F, na.rm = T)
    
    line_df.temp$scales = sapply(1:dim(line_df.temp)[1], FUN = function(i){
      var.temp = as.character( line_df.temp$variable )
      line_df.temp[i, paste0('scales', substr(var.temp[i], nchar(var.temp[i])-nlen+1, nchar(var.temp[i])))]
    })
    
    line_df = line_df.temp[, (numObj+1):(numObj+5)]
    
    # 2. data frame that is used to plot CI
    o1 = obj[,c(scales_names, 'dataset_h', 'dataset_v', low_names)] 
    o2 = obj[,c(scales_names, 'dataset_h', 'dataset_v', high_names)] 
    
    low_df = melt(o1, measure.vars = low_names, variable.name = 'low', factorsAsStrings = F, na.rm = T )
    high_df = melt(o2, measure.vars = high_names, variable.name = 'high', factorsAsStrings = F, na.rm = T )
    
    colnames(low_df) = c(scales_names, 'dataset_h', 'dataset_v', 'low', 'l_value')
    CI_df.temp =  data.frame(low_df, 
                             high = high_df$high, 
                             h_value = high_df$value) #levels(CI_df$dataset_h), levels(CI_df$dataset_v) no need to set again
    
    CI_df.temp$scales = sapply(1:dim(CI_df.temp)[1], FUN = function(i){
      var.temp = as.character( CI_df.temp$low )
      CI_df.temp[i, paste0('scales', substr(var.temp[i], nchar(var.temp[i])-nlen+1, nchar(var.temp)))]
    })
    CI_df = CI_df.temp[, (numObj+1):(numObj+7)]
    
    
    # S4: Generate the graph
    # CALL Graphical Functions
    p = ggplot() + 
      geom_line( data = line_df, mapping = aes(x = scales, y = value, color = variable, linetype = variable) ) + 
      geom_point(data = line_df, mapping = aes(x = scales, y = value, color = variable, size = variable, shape = variable) ) +
      
      geom_ribbon(data = CI_df, mapping = aes(x = scales, ymin = l_value, ymax = h_value, fill = low),alpha = transparence ) +
      
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) + 
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
      
      scale_linetype_manual(values = c(line.type)) +
      scale_shape_manual(values = c(point.shape))+
      scale_size_manual(values = c(point.size)) +
      scale_color_manual(values = c(line.color)) + 
      scale_fill_manual(values = c(CI.color)) 
    
    if( background == 'white' ){
      p = p + theme_bw() 
    }
    
    #p = p + facet_grid(dataset_h~dataset_v, labeller = label_parsed) + 
    p = p + facet_grid(dataset_h~dataset_v, labeller = label_parsed) + 
      theme(legend.position='none') + 
      theme(strip.background = element_rect(fill= facet.label.background) )
    
    p = p +
      xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
      theme(
        plot.title = element_text(size= title.size),
        axis.title.y = element_text(size= axis.label.size),
        axis.text.y  = element_text(size= axis.tick.size),
        axis.title.x = element_text(size= axis.label.size),
        axis.text.x  = element_text(size= axis.tick.size),
        strip.text = element_text(size = facet.label.size) )
    
    p
    
  } # end else
}

#' Graphically Compare Classical with Robust GMWM Models
#'
#' Creates a table of graphs to compare GMWM models constructed by classical and robust methods.
#' @param ...          Several \code{gmwm} objects, and they must be constrcuted by the same data and same model. Or one \code{rgmwm} object created by \code{\link{rgmwm}} function.
#' @param display.eff  A \code{boolean} indicating whether the classical/robust method with efficiency should be displayed in the facet label. It is only used when \code{facet.label} is NULL.
#' @param order.eff    A \code{boolean} indicating whether the object should be ordered by their efficiency value.
#' @param facet.label  A \code{character vector} indicating what should be displayed in the facet labels.
#' @param background   A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} between 0 to 1 that controls the transparency of confidence interval.
#' @param CI.color     A \code{character vector} that indicates the color of the confidence interval (e.g. 'black', 'red', '#003C7D', etc.)
#' @param line.type    A \code{character vector} that indicates the type of lines.
#' @param line.color   A \code{character vector} that indicates the color of lines.
#' @param point.size   A \code{integer vector} that indicates the size of points on lines. 
#' @param point.shape  A \code{integer vector} that indicates the shape of points on lines.
#' @param title        A \code{string} that indicates the title of the graph.
#' @param title.size   An \code{integer} that indicates the size of title.
#' @param axis.label.size        An \code{integer} that indicates the size of label.
#' @param axis.tick.size         An \code{integer} that indicates the size of tick mark.
#' @param facet.label.size       An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @author Wenchao
#' 
#' @section Constraints on \code{gmwm} objects:
#' This function is designed to compare the \code{gmwm} objects with different computation methods and efficiency values. Therefore,
#' \code{gmwm} objects supplied to this function are assumed to have different efficiency values and are constructed by 
#' the same data and same model. The function will check the constraint before generating the graph.
#' 
#' @section How to change \code{CI.color} and \code{facet.label}:
#' To change \code{CI.color} and \code{facet.label}, the user must supply a vector. 
#' If the user compares \code{N} {gmwm} objects, then he/she must supply a vector of length 
#' \code{N} to \code{CI.color} and \code{facet.label}. Please check examples for help. Additionally,
#' \code{order.eff} will change the order how objects are plotted. User is recommended to plot the graph
#' without changing the aesthetics firstly. Then, the user can generate the graph by supplying elements to graphical parameters, e.g. \code{CI.color},
#' \code{facet.label}, in left-to-right order.
#' 
#' @section How to change other graphical aesthetics:
#' \code{line.type}, \code{line.color}, \code{point.size}, \code{point.size} must be a \code{vector}. 
#' If the user compares \code{N} \code{gmwm} objects, then a vector of length \eqn{4 \times N}{ 4 x N } must be supplied to these 
#' parameters. For example, if the user wants to compare \code{gmwm1} with \code{gmwm2}, the order to change the 
#' graphical aesthetics is:
#' "gmwm1 Empirical WV, gmwm2 Empirical WV, gmwm1 lower bound of CI, gmwm2 lower bound of CI, gmwm1 higher bound of 
#' CI, gmwm2 higher bound of CI, gmwm1 Implied WV, gmwm2 Implied WV."
#' 
#' Please check examples for help.
#' 
#' 
#' @section Other details:
#' 
#' \code{melt()} is used to convert the data frame from wide format to long format. However,
#' \code{melt()} cannot deal with expressions, so the user cannot supply expressions to the parameter 
#' \code{facet.label}. Users should convert expressions into characters via \code{as.character()}.
#' 
#' If you meet the error "polygon edge not found", it is complaining that you don't have enough space to
#' plot the graph. If you are using RStudio, you can adjust the plot window. You can also open graphical window by using
#' \code{quartz()} (Mac/Linux) or \code{windows()} (Windows PC).
#' 
#' @examples
#' \dontrun{
#' #Simulate data
#' set.seed(8836)
#' n = 1000
#' x = gen_gts(n, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' x[1] = 1000 #Robust gmwm works better for contaminated data
#' GMWM1 = gmwm(2*AR1()+RW(), data = x, robust  = TRUE, eff = 0.9)
#' GMWM2 = gmwm(2*AR1()+RW(), data = x, robust  = TRUE, eff = 0.6)
#' GMWM3 = gmwm(2*AR1()+RW(), data = x, robust  = FALSE)
#' 
#' # How order.eff works:
#' compare_eff(GMWM1, GMWM2, GMWM3, order.eff = FALSE)
#' compare_eff(GMWM1, GMWM2, GMWM3, order.eff = TRUE) 
#' # order.eff=TRUE: The objects are plotted in descending order of model efficiency.
#'
#' 
#' # How to modify facet.label
#' compare_eff(GMWM2, GMWM3, order.eff = FALSE, facet.label = c('Rob. eff. 0.6', 'Cla.') )
#' compare_eff(GMWM2, GMWM3, order.eff = TRUE, facet.label = c('Cla.', 'Rob. eff. 0.6') )
#' 
#' # How to modify graph aesthetics
#' compare_eff(GMWM2, GMWM3, order.eff = FALSE, CI.color = c("#003C7D", "#F47F24"),
#'            line.color = rep(c("#003C7D", "#F47F24"), 4) )
#' compare_eff(GMWM2, GMWM3, order.eff = TRUE, CI.color = c("#F47F24", "#003C7D"),
#'            line.color = rep(c("#F47F24", "#003C7D"), 4) )
#' }
compare_eff = function(..., display.eff = T, order.eff = T, facet.label = NULL, 
                       background = 'white', transparence = 0.05, CI.color = NULL,
                       line.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                       title = NULL, title.size= 18, 
                       axis.label.size = 16, axis.tick.size = 11, 
                       facet.label.size = 13, facet.label.background = "#003C7D33",
                       axis.x.label = expression(paste("Scale ", tau)),
                       axis.y.label = expression(paste("Wavelet Variance ", nu))){
  scales=variable=l_value=h_value=low=.x=NULL
  
  # S1: Checking statement (Reset it to default setting if user passes wrong values)
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...) 
  numObj = length(obj_list)
  object.names = as.character(substitute(...()))
  
  #deal with rgmwm object
  if( is(obj_list[[1]], 'rgmwm') ){
    if(numObj == 1){
      obj_list = obj_list[[1]]
      numObj = length(obj_list)
    }else{
      stop("Please supply only one 'rgmwm' object.")
    }
  }
  
  if(numObj == 1){
    #plot.gmwm will do the parameter checking
    #other parameters are not listed here, and they cannot be passed to plot.gmww by '...'
    warning('One object is supplied. You are actually calling plot.gmwm().')
    
    if(is.null(CI.color)) {CI.color="#003C7D"}
    
    plot.gmwm(x = obj_list[[1]], process.decomp = F, background = background, transparence = transparence, CI.color = CI.color, CI = T, bw = F, line.type = line.type,
              line.color = line.color, point.size = point.size, point.shape = point.shape, title = title,
              title.size = title.size, axis.label.size = axis.label.size, axis.tick.size = axis.tick.size,
              axis.x.label = axis.x.label, axis.y.label = axis.y.label)
    
  }else{
    # check whter the objects supplied is valid
    # if not, stop
    is.validCompEffObj(obj_list)
    
    # Re-order the object list
    if(order.eff){
      eff.vec = getEff(obj_list)
      
      # if it is in descending order, it's ok
      if(is.unsorted(-eff.vec)){
        new.order = order(eff.vec, decreasing = T)
        
        obj_list = obj_list[new.order]
        object.names = object.names[new.order]
      }
    }
    
    #check parameter
    params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape', 'facet.label')
    requireLength = c(4*numObj, 4*numObj, numObj, 4*numObj, 4*numObj, numObj)
    default = list(NULL, NULL, NULL, NULL, NULL, NULL)
    nullIsFine = c(rep(T,6))

    checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
    
    # S2: Auto-select parameters, if not provided by users
    
    # decide what should appear in facet label
    if(is.null(facet.label) && display.eff){
      object.names = sapply(obj_list, FUN = function(x){formatRobustEff(x)} )
    }else if(!is.null(facet.label)){
      object.names = facet.label 
    }
    
    # Problem: When object names are the same, function will not work
    # Solution: Add space after object names
    object.names = addSpaceIfDuplicate(object.names)
    
    # title
    if(is.null(title)){
      
      robustMod = sapply(obj_list, FUN = function(x){x$robust})
      if(any(robustMod == F)){ #classical model exists
        title = 'Comparison of Classical vs. Robust GMWM Models'
      }else{
        title = 'Comparison of Robust GMWM Models'
      }
    }
    
    # in the order: WV low high theo
    if( is.null(line.type) ){line.type = c(rep('solid',numObj),#WV
                                         rep('dotted', numObj),#low
                                         rep('dotted', numObj),#high
                                         rep('solid', numObj))}#theo
    
    if(is.null(CI.color)){
      CI.color = switch(as.character(numObj),
                        '2' = c("#003C7D","#F47F24"), #UIUC
                        '3' = c("#003C7D","#F47F24","#34A853"), 
                        ggColor(numObj) )
                        
    }
    
    Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
    if(is.null(line.color)) {
      if (numObj <= 3) {
        line.color = switch(as.character(numObj),
                            '2' = c(rep(CI.color,3), c("#E41A1C","#4DAF4A")),#set 1
                            '3' = c(rep(CI.color,3), c("#E41A1C", "#4DAF4A", "#341456"))) #set 1 and purple
      }else{
        modulus = numObj %/% 8
        remainder = numObj %% 8
        theo.color =  c(rep(Set1, times = modulus), Set1[1:remainder])
        line.color = c(rep(CI.color, 3), #WV,low,high
                       theo.color) #theo
        if(numObj>8){
          warning('More than 8 objects are supplied, but the palette has only 8 colors')
        }
       
      }
    }
    
    #point.size will auto-decrease as the number of models increases
    if(is.null(point.size)){
      
      if(numObj > 10){temp.point.size = 1
      }else{
        temp.point.size = 4.5-0.4*(numObj-1)}
      
      point.size = c(rep(temp.point.size, numObj),#WV
                     rep(0, numObj),#low
                     rep(0, numObj),#high
                     rep(temp.point.size, numObj))#theo
      
    }
    
    if(is.null(point.shape)){point.shape = c(rep(20, numObj),#WV
                                             rep(46, numObj),#low
                                             rep(46, numObj),#high
                                             rep(1, numObj))}#theo
    
    # S3: Rearrange the data into a data frame which can be passed to next step
    total.len = 0
    each.len = numeric(numObj)
    
    for (i in 1:numObj ){
        # freq conversion
        obj_list[[i]]$scales = obj_list[[i]]$scales/obj_list[[i]]$freq
        
        each.len[i] = length(obj_list[[i]]$wv.empir)
        total.len = total.len + each.len[i]*numObj
    }
    
    #Choose the column names for empty data frame
    nlen = nchar(numObj)
    WV_names = sprintf(paste0("WV%0",nlen,"d"),1:numObj)
    low_names = sprintf(paste0("low%0",nlen,"d"),1:numObj)
    high_names = sprintf(paste0("high%0",nlen,"d"),1:numObj)
    z_theo_names = sprintf(paste0("z_theo%0",nlen,"d"),1:numObj)
    col.names = c('scales', 'dataset_h', 'dataset_v', 
                  WV_names, low_names, high_names, z_theo_names)
    
    #Initialize empty data frame with right number of rows
    num.col = 3 + 4*numObj
    obj = as.data.frame(matrix(NA, ncol = num.col, nrow = total.len), stringsAsFactors = FALSE)
    colnames(obj) = col.names
    
    #put data into data frame
    t = 1
    for (i in 1:numObj){
      for(j in 1:numObj){
        
        ## Diagonal ========================
        if(i == j){
          #compare to itself
          #set some values to NA
          d = each.len[i]
          
          WV_name = WV_names[i]
          low_name = low_names[i] 
          high_name = high_names[i] 
          z_theo_name = z_theo_names[i]
          
          obj[t:(t+d-1),c('scales', 'dataset_h', 'dataset_v', 
                          WV_name, low_name, high_name, z_theo_name)] = 
            
            data.frame(obj_list[[i]]$scales,
                       
                       object.names[i],
                       object.names[j],
                       
                       obj_list[[i]]$wv.empir,
                       obj_list[[i]]$ci.low,
                       obj_list[[i]]$ci.high,
                       obj_list[[i]]$theo, stringsAsFactors = F)
          
          t = t +d
          
        }else if (i > j ){
          ## lower triagular ========================
          
          d = each.len[i]
          
          WV_name = WV_names[c(i,j)]
          low_name = low_names[c(i,j)] 
          high_name = high_names[c(i,j)] 
          #z_theo_name = z_theo_names[c(i,j)]
          
          obj[t:(t+d-1),c('scales', 'dataset_h', 'dataset_v', 
                          WV_name, low_name, high_name)] = 
            
            data.frame(obj_list[[i]]$scales,
                       
                       object.names[i],
                       object.names[j],
                       
                       obj_list[[i]]$wv.empir,
                       obj_list[[j]]$wv.empir,
                       
                       obj_list[[i]]$ci.low,
                       obj_list[[j]]$ci.low,
                       
                       obj_list[[i]]$ci.high,
                       obj_list[[j]]$ci.high, stringsAsFactors = F)
          
          t = t +d
          
        }else{
          ## upper triagular ======================== 
          
          d = each.len[i]
          
          z_theo_name = z_theo_names[c(i,j)]
          
          obj[t:(t+d-1),c('scales','dataset_h','dataset_v',z_theo_name)] = 
            
            data.frame(obj_list[[i]]$scales,
                       
                       object.names[i],
                       object.names[j],
                       
                       obj_list[[i]]$theo,
                       obj_list[[j]]$theo, stringsAsFactors = F)
          
          t = t +d
        }
      } # end of j loop
    }# end of i loop
    
    # change the order of facet label
    obj$dataset_h = factor(obj$dataset_h, levels = object.names )
    obj$dataset_v = factor(obj$dataset_v, levels = object.names )
    
    # 1. data frame that is used to plot lines
    line_df = melt(obj, id.vars = c('scales', 'dataset_v', 'dataset_h'), factorsAsStrings = F, na.rm = T)
    
    # 2. data frame that is used to plot CI
    o1 = obj[,c('scales', 'dataset_h', 'dataset_v', low_names)] 
    o2 = obj[,c('scales', 'dataset_h', 'dataset_v', high_names)] 
    
    low_df = melt(o1, measure.vars = low_names, variable.name = 'low', factorsAsStrings = F, na.rm = T )
    high_df = melt(o2, measure.vars = high_names, variable.name = 'high', factorsAsStrings = F, na.rm = T )
    
    colnames(low_df) = c('scales', 'dataset_h', 'dataset_v', 'low', 'l_value')
    CI_df =  data.frame(low_df, 
                        high = high_df$high, 
                        h_value = high_df$value) #levels(CI_df$dataset_h), levels(CI_df$dataset_v) no need to set again
    
    # S4: Generate the graph
    # CALL Graphical Functions
    p = ggplot() + 
      geom_line( data = line_df, mapping = aes(x = scales, y = value, color = variable, linetype = variable) ) + 
      geom_point(data = line_df, mapping = aes(x = scales, y = value, color = variable, size = variable, shape = variable) ) +
      
      geom_ribbon(data = CI_df, mapping = aes(x = scales, ymin = l_value, ymax = h_value, fill = low),alpha = transparence ) +
      
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) + 
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
    
      scale_linetype_manual(values = c(line.type)) +
      scale_shape_manual(values = c(point.shape))+
      scale_size_manual(values = c(point.size)) +
      scale_color_manual(values = c(line.color)) + 
      scale_fill_manual(values = c(CI.color)) 
    
    if( background == 'white' ){
      p = p + theme_bw() 
    }
    
    #p = p + facet_grid(dataset_h~dataset_v, labeller = label_parsed) + 
    p = p + facet_grid(dataset_h~dataset_v) + 
      theme(legend.position='none') + 
      theme(strip.background = element_rect(fill= facet.label.background) )
    
    p = p +
      xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
      theme(
        plot.title = element_text(size= title.size),
        axis.title.y = element_text(size= axis.label.size),
        axis.text.y  = element_text(size= axis.tick.size),
        axis.title.x = element_text(size= axis.label.size),
        axis.text.x  = element_text(size= axis.tick.size),
        strip.text = element_text(size = facet.label.size) )
   
    p
    
  } # end else
}

