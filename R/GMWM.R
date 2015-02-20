

#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param desc A \code{vec<string>} containing one of the allowed models. Not used if \code{type="arma"}
#' @param wvcov A \code{wvcov} object
#' @param signal The time series being studied
#' @param type A \code{string} containing the type of GMWM needed e.g. IMU, ARMA, or SSM
#' @param params A \code{list} containing the numbers of AR and MA parameters. Only used if \code{type="arma"}
#' @param B A \code{integer} to sample the space for IMU and SSM models to ensure AR1 identitability.
#' @name gmwm_guide
#' @examples
#' # ARMA case
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sig2 = 1) + gen_ar1(n,phi=0.95, sig2 = .1)
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save = gmwm(desc=c("AR1"), wvcov=out, signal=x, type="imu", B = 10000)
#'  
#' # ARMA case
#' set.seed(1336)
#' x=arima.sim(n = 200, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796))
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
#' save = gmwm(wvcov=out, signal=x, type="arma", params=list(ar=2,ma=2))
gmwm = function(desc, wvcov, signal, type="imu", params = list(ar = 1, ma = 2), B = 1000){
  
  if(type == "imu" || type == "ssm"){
    if(length(desc) == 0 ){
      stop("desc must contain a list of components (e.g. AR1, DR, RW, QN, WN). Do NOT use params!")
    }
    
    np = num_model_params(desc)
    
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }
    
    if(wvcov$robust){
      np = np+1
      if(np > length(wvcov$scales)){
        stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
      }
      ind = 1:np
      wvcov$scales = wvcov$scales[ind]
      wvcov$wv.empir = wvcov$wv.empir[ind]
      wvcov$V = wvcov$V[ind,ind]
      wvcov$ci_low = wvcov$ci_low[ind]
      wvcov$ci_high = wvcov$ci_high[ind]
    }
    
    N = length(signal)
    
    if(type=="imu"){
      out = .Call('GMWM_gmwm_imu_cpp', PACKAGE = 'GMWM', desc, signal, wvcov$V, wvcov$wv.empir, wvcov$scales, N, B)
    }else{
      out = .Call('GMWM_gmwm_ssm_cpp', PACKAGE = 'GMWM', desc, signal, wvcov$V, wvcov$wv.empir, wvcov$scales, N, B)
    }
    theo = theo_wv(out, desc, wvcov$wv.empir, wvcov$scales, N)
  } else if(type == "arma"){
    if(length(params) != 2){
      stop("params must be a list with a length of two. The AR component is first and the MA component is second.")
    }
    p = params[[1]]
    q = params[[2]]
    np = p+q+1
    if(np > length(wvcov$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
    }

    if(wvcov$robust){
      np = np+1
      if(np > length(wvcov$scales)){
        stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need the same number of scales as parameters to estimate.")
      }
      ind = 1:np
      wvcov$scales = wvcov$scales[ind]
      wvcov$wv.empir = wvcov$wv.empir[ind]
      wvcov$V = wvcov$V[ind,ind]
      wvcov$ci_low = wvcov$ci_low[ind]
      wvcov$ci_high = wvcov$ci_high[ind]
    }
    
    guess.start = arima(signal, order=c(p,0,q), include.mean=FALSE, method="CSS")
    theta = c(coef(guess.start), guess.start$sigma2)
    
    out = .Call('GMWM_gmwm_arma_cpp', PACKAGE = 'GMWM', theta, wvcov$V, p,  q, wvcov$scales, wvcov$wv.empir)
    

    if(p == 0){
      ar = numeric()
    }else{
      ar = out[1:p]
    }
    
    if(q == 0){
      ma = numeric()
    }else{
      ma = out[(p+1):(length(out)-1)]
    }

    theo = arma_to_wv(ar, ma, wvcov$scales, out[length(out)])
  } 
    
  out = structure(list(estimate = as.vector(out),
                       wv.empir=wvcov$wv.empir, 
                       ci_low=wvcov$ci_low, 
                       ci_high = wvcov$ci_high, 
                       scales = wvcov$scales, 
                       theo = theo, 
                       robust = wvcov$robust,
                       eff = wvcov$eff,
                       type = type), class = "gmwm")
  invisible(out)
}


adv.gmwm = function(desc, theta, wvcov){
  .Call('GMWM_adv_gmwm_cpp', PACKAGE = 'GMWM', theta, desc, V, wv_empir, tau, N)
}


autoplot.gmwm = function(x){
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = x$wv.empir, low = x$ci_low, high = x$ci_high, scale = x$scales, theo = x$theo)
  CI = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3")) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Haar Wavelet Variance Representation") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta),")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0)) 
  
  CI
}

autoplot.comp = function(x){
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#EF8A1C","LINE4"="#6E8BBF")
  
  WV = data.frame(emp = x$wv.empir, low = x$ci_low, high = x$ci_high, scale = x$scales, theo = x$theo1, model2 = x$theo2)
  CI = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3"), size = 1) + 
    geom_line(aes(y = model2, colour = "LINE4"), size = 1) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    geom_point(aes(y = model2, colour = "LINE4"), size = 4, shape = 0) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Haar Wavelet Variance Representation") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Model 1: ", nu,"(",hat(theta)[1],")")),expression(paste("Model 2: ", nu,"(",hat(theta)[2],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1,1), colour = c("#000000","#999999","#EF8A1C","#6E8BBF")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0)) 
  
  CI
}



autoplot.comp2 = function(x){
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = x$wv.empir, low = x$ci_low, high = x$ci_high, scale = x$scales, theo = x$theo1, model2 = x$theo2)
  
  CI1 = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3")) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Model 1") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta)[1],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0))
  
  CI2 = ggplot(WV, aes( x = scale, y = low)) + geom_line(aes(colour = "LINE2"), linetype = "dotted") +
    geom_line(aes(y = high, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = emp, colour = "LINE1")) + geom_point(aes(y = emp, colour = "LINE1"), size = 3) +
    geom_line(aes(y = model2, colour = "LINE3")) + 
    geom_point(aes(y = model2, colour = "LINE3"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) +
    ggtitle("Model 2") + 
    scale_colour_manual(name=" ", labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                           expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta)[2],")"))), 
                        values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
    guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9")))) +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0))
  
  multiplot(CI1, CI2, cols=2)	
}

compare.models = function(GMWM1, GMWM2, split = FALSE){
  x = data.frame(wv.empir = GMWM1$wv.empir, ci_low = GMWM1$ci_low, 
                 ci_high = GMWM1$ci_high, scales = GMWM1$scales, theo1 = GMWM1$theo, theo2 = GMWM2$theo) 
  if (split == TRUE){
    class(x) = "comp2"
  }else{
    class(x) = "comp"
  }
  autoplot(x)
}
