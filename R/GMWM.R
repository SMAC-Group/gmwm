#' @title GMWM for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @param data A \code{matrix} or \code{data.frame} object with only column (e.g. \eqn{N \times 1}{ N x 1 })
#' @param model.type A \code{string} containing the type of GMWM needed e.g. IMU or SSM
#' @param compute.v A \code{string} indicating the type of covariance matrix solver. "fast", "bootstrap", "asymp.diag", "asymp.comp", "fft"
#' @param augmented A \code{boolean} indicating whether to add additional moments (e.g. mean for drift and variance for all other components).
#' @param alpha A \code{double} between 0 and 1 that correspondings to the \eqn{\frac{\alpha}{2}}{alpha/2} value for the wavelet confidence intervals.
#' @param robust A \code{boolean} indicating whether to use the robust computation (TRUE) or not (FALSE).
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency.
#' @param G An \code{integer} to sample the space for IMU and SSM models to ensure optimal identitability.
#' @param K An \code{integer} that controls how many times the bootstrapping procedure will be initiated.
#' @param H An \code{integer} that indicates how many different samples the bootstrap will be collect.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' What is NOT active:
#' 1. Augmented GMWM (additional moments)
#' 2. Guessing ARMA models is NOT supported. You must supply specific parameters (see example below)
#' 3. ARMA Bootstrapping.
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
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' data = gen.ts(AR1(phi = .99, sigma2 = 0.01) + WN(sigma2 = 1), n)
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
#' data = arima.sim(n = 200, 
#'               list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
#'               sd = sqrt(0.1796))
#' #guided.arma = gmwm(ARMA(2,2), data, model.type="ssm")
#' adv.arma = gmwm(ARMA(ar=c(0.8897, -0.4858), ma = c(-0.2279, 0.2488), sigma2=0.1796),
#'                 data, model.type="ssm")
gmwm = function(model, data, model.type="ssm", compute.v="auto", augmented=FALSE, robust=FALSE, eff=0.6, inference = NULL, model.select = NULL, alpha = 0.05, seed = NULL, G = NULL, K = 1, H = 100){
  
  # Are we receiving one column of data?
  if( (class(data) == "data.frame" && ncol(data) > 1) || ( class(data) == "matrix" && ncol(data) > 1 ) ){
    stop("The function can only operate on one column of data")
  }
  
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), AR(p), MA(q), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  np = model$plength
  
  N = length(data)
  
  starting = model$starting
  
  # Input guessing
  if((is.null(G) & starting) || !is.whole(G)){
    if(N > 10000){
      G = 1000
    }else{
      G = 20000
    }
  }else if(!starting){
    G = 0
  }
  
  # For reproducibility
  if(is.null(seed) || !is.whole(seed)){
    seed = floor(runif(1, 1, 10000))
  }
  
  set.seed(seed)
  
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  # Model type issues
  if(model.type != "imu" && model.type != "ssm"){
    stop("Model Type must be either IMU or SSM!")
  }
  
  # Verify Scales and Parameter Space
  nlevels =  floor(log2(length(data)))
  scales = .Call('GMWM_scales_cpp', PACKAGE = 'GMWM', nlevels)
  
  if(np > length(scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need at least the same number of scales as parameters to estimate.")
  }
  
  if(robust){
    np = np+1
    if(np > length(scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we  at least need the same number of scales as parameters to estimate.")
    }
  }
  
  # Compute inference on small time series.
  if(is.null(inference)){
    if(N > 10000){
      inference = FALSE
    }else{
      inference = TRUE
    } 
  }
  
  
  # Compute model.select on small time series.
  if(is.null(model.select)){
    if(N > 10000){
      model.select = FALSE
    }else{
      model.select = TRUE
    } 
  }
  
  # Auto setting

  # Compute fast covariance if large sample, otherwise, bootstrap.
  if(compute.v == "auto" || ( compute.v != "fast" && compute.v != "bootstrap")){
    if(N > 10000){
      compute.v = "fast"
    }else{
      compute.v = "bootstrap"
    }
  }
  

  theta = model$theta

  out = .Call('GMWM_gmwm_master_cpp', PACKAGE = 'GMWM', data, theta, desc, obj, model.type, starting = model$starting,
                                                         p = alpha, compute_v = compute.v, K = K, H = H, G = G,
                                                         robust=robust, eff = eff, inference || model.select)
  #colnames(out) = model$desc
  
  estimate = out[[1]]
  rownames(estimate) = model$process.desc
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  
  out = structure(list(data = data,
                       estimate = estimate,
                       init.guess = init.guess,
                       wv.empir = out[[3]], 
                       ci.low = out[[4]], 
                       ci.high = out[[5]],
                       V = out[[6]],
                       omega = out[[12]],
                       obj.fun = out[[11]],
                       orgV = out[[7]],
                       theo = out[[9]],
                       decomp.theo = out[[10]],
                       scales = scales, 
                       robust = robust,
                       eff = eff,
                       model.type = model.type,
                       compute.v = compute.v,
                       augmented = augmented,
                       alpha = alpha,
                       expect.diff = out[[8]],
                       N = N,
                       G = G,
                       H = H,
                       K = K,
                       model = model,
                       starting = model$starting,
                       inference = inference,
                       model.select = model.select,
                       seed = seed), class = "gmwm")
  invisible(out)
}

#' @title Update GMWM object for IMU, ARMA, SSM, and Robust
#' @description GMM object
#' @param object A \code{gmwm} object.
#' @param model A \code{ts.model} object containing one of the allowed models.
#' @return A \code{gmwm} object that contains:
#' \itemize{
#'  \item{}
#'  \item{}
#'  \item{}
#' }
#' @details
#' This function is under work. Some of the features are active. Others... Not so much. 
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' exact.model = AR1(phi=.99, sigma2 = 0.01) + WN(sigma2=1)
#' data = gen.ts(exact.model)
#' 
#' # Create an initial model that is not accurate
#' bad.model = gmwm(AR1(), data = data)
#' 
#' # Models can contain specific parameters e.g.
#' updated.model = update(bad.model, exact.model)
#' 
#' # Or...
#' updated.model.guided = update(bad.model, AR1()+AR1())
update.gmwm = function(object, model, ...){
  # Do we have a valid model?
  if(!is(model, "ts.model")){
    stop("model must be created from a ts.model object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  np = model$plength
  
  # Information used in summary.gmwm:
  summary.desc = model$desc
  
  models.active = count_models(desc)
  # Identifiability issues
  if(any( models.active[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  # ID error:
  if( sum(models.active) == 1 & models.active["ARMA"] == 1 & model$starting){
    warning("ARMA starting guesses using update.gmwm are NOT based on CSS but an alternative algorithm.")
  }
  
  if(np > length(object$scales)){
    stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need  at least  the same number of scales as parameters to estimate.")
  }
  
  if(object$robust){
    np = np+1
    if(np > length(object$scales)){
      stop("Please supply a longer signal / time series in order to use the GMWM. This is because we need one additional scale since robust requires the amount of parameters + 1 to estimate.")
    }
  }
  
  out = .Call('GMWM_gmwm_update_cpp', PACKAGE = 'GMWM',
                  model$theta,
                  desc, obj, 
                  object$model.type, object$N, object$expect.diff, 
                  object$orgV, object$scales, object$wv.empir,
                  model$starting, 
                  object$compute.v, object$K, object$H,
                  object$G, 
                  object$robust, object$eff, object$inference || object$model.select)

  estimate = out[[1]]
  rownames(estimate) = model$process.desc
  init.guess = out[[2]]
  rownames(init.guess) = model$process.desc
  
  object$estimate = estimate
  object$init.guess = init.guess
  
  object$V = out[[3]]
  object$theo = out[[4]]
  object$decomp.theo = out[[5]]
  
  object$starting = model$starting

  invisible(object)
}

#' @title Summary of GMWM object
#' @description Displays summary information about GMWM object
#' @method summary gmwm
#' @param object A \code{GMWM} object
#' @param ci.bootstrap A value containing either: NULL (auto), TRUE, or FALSE
#' @param model.select A value containing either: NULL (auto), TRUE, FALSE
#' @param ... other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen.ts(AR1(phi=.1, sigma2 = 1) + AR2(phi=0.95, sigma2 = .1),n)
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' summary(mod)
summary.gmwm = function(object, ci.bootstrap = NULL, B = 20, ...){
  out = cbind(object$init.guess, object$estimate)
  colnames(out) = c("Initial Guess", "Estimates") 
  
  inference = object$inference
  
  model.select = object$model.select

  if(is.null(ci.bootstrap)){
    if(object$N > 10000){
      ci.bootstrap = FALSE
    }else{
      ci.bootstrap = TRUE
    } 
  }
  
  if("ARMA" %in% object$model$desc){
    
    if( model.select == T){
      warning("ARMA is not currently supported for model selection. The model selection results will not be displayed.")
      model.select = FALSE
    }
    
    if(ci.bootstrap == FALSE){
      warning(paste0("The numerical derivative of ARMA(p,q), where p > 1 and q > 1, may be inaccurate leading to inappropriate CIs.\n",
              "Consider using the ci.bootstrap = T option on the summary function."))
    }
  }
  
  
  if(inference || model.select){

    mm = .Call('GMWM_get_summary', PACKAGE = 'GMWM', object$estimate, object$model$desc, object$model$obj.desc,
                                              object$model.type, object$wv.empir, object$theo, 
                                              object$scales, object$V, object$omega, object$N, object$alpha, inference, model.select, ci.bootstrap, B, object$robust, object$eff)
  }else{
    mm = vector('list',3)
    mm[1:3] = NA
  }
  
  if(inference){
    out.coln = colnames(out)
    out = cbind(out, mm[[1]])
    colnames(out) = c(out.coln, "CI Low", "CI High", "SE")
  }
  
  x = structure(list(estimates=out, 
                     testinfo=mm[[2]],
                     model.score = mm[[3]],
                     inference = inference, 
                     model.select = model.select, 
                     starting = object$starting,
                     seed = object$seed,
                     N = object$N), class = "summary.gmwm")
}

#' @title Print summary.gmwm object
#' @description Displays summary information about GMWM object
#' @method print summary.gmwm
#' @param x A \code{GMWM} object
#' @param ... other arguments passed to specific methods
#' @return Text output via print
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' xt = gen.ts(AR1(phi=.1, sigma2 = 1) + AR2(phi=0.95, sigma2 = .1),n)
#' mod = gmwm(AR1()+AR1(), data=xt, model.type="imu")
#' summary(mod)
print.summary.gmwm = function(x, ...){
  cat("Model Information: \n")
  print(x$estimates)
  if(x$starting){
    cat("\nThe values supplied under initial guess were generated by the program.\n")
    cat(paste0("To replicate the guessed parameters, use seed: ",x$seed, "\n"))
  }else{
    cat("\nThe values supplied under initial guess were given by YOU!\n")
  }
  
  if(x$inference){
    cat(paste0("\nThe Goodness of Fit test statistic is: ", round(x$testinfo[1],2),
          " on ",x$testinfo[3]," degrees of freedom\n",
          "The resulting p-value is: ", round(x$testinfo[2],4)),"\n\n")
  }else{
    cat("\nInference was not run. \nTo obtain theta confidence intervals and Goodness of Fit information, please re-run summary() with inference = TRUE.")
  }
  
  if(x$model.select){
    cat(paste0("\nThe model score statistic is: ", round(x$model.score[[1]]),"\n\n"))
  }
}


#' @title Predict future points in the time series using the solution of the Generalized Method of Wavelet Moments
#' @description Creates a prediction using the estimated values of GMWM through the ARIMA function within R.
#' @method predict gmwm
#' @param object GMWM 
#' @param n.ahead Number of observations to guess.
#' @param data Set used to create the GMWM.
predict.gmwm = function(object, n.ahead, ...){
  
  ts.mod = object$model
  
  if(length(ts.mod$desc) > 1 || ts.mod$desc != "ARMA")
    stop("The predict function only works with stand-alone ARMA models.")
  
  objdesc = ts.mod$obj.desc[[1]]
  
  p = objdesc[1]
  q = objdesc[2]
  
  mod = arima(object$data, order = c(p, 0, q),
              method="ML",
              fixed = object$estimate[1:(p+q)],
              transform.pars = F,
              include.mean = F)
  
  pred = predict(mod, n.ahead = n.ahead, newxreg = NULL,
                 se.fit = TRUE, ...)
  
  
  out = structure(list(pred = pred$pred,
                       se = pred$se,
                       resid = mod$residuals)
                  , class = "predict.gmwm")
                          
}

#' @title Wrapper to Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method plot gmwm
#' @param x A \code{GMWM} object
#' @param individual A \code{boolean} that indicates whether the graph should be plotted individually or not
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param color.CI A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param ... other arguments passed to specific methods
#' @note Use \code{?autoplot.gmwm} and \code{?autoplot.gmwm2} to view specific explanation for parameters \code{color.line}, \code{line.type}. They are different when plotting individually or not 
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen_ar1(n, phi=.1, sigma2 = 1) + gen_ar1(n,phi=0.95, sigma2 = .1)
#' mod = gmwm(AR1(), data=x, model.type="imu")
#' plot(mod)
plot.gmwm = function(x, individual = FALSE, CI = T, transparence = 0.1, 
                     color.CI = "#003C7D",  
                     graph.title = NA, graph.title.size= 15, 
                     axis.label.size = 13, axis.tick.size = 11, 
                     title.x.axis = expression(paste("Scale ", tau)),
                     title.y.axis = expression(paste("Wavelet Variance ", nu)),
                     legend.title = '', legend.key.size = 1, legend.title.size = 13, 
                     legend.text.size = 13, ...){
  
  if(individual){
    # plot individually
    autoplot.gmwm2(x, CI = CI, transparence = transparence, 
                   color.CI = color.CI,  
                   graph.title = graph.title, graph.title.size= graph.title.size, 
                   axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                   title.x.axis = title.x.axis,
                   title.y.axis = title.y.axis,
                   legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                   legend.text.size = legend.text.size)
  }
  else{
    autoplot.gmwm(x, CI = CI, transparence = transparence, 
                  color.CI = color.CI,  
                  graph.title = graph.title, graph.title.size= graph.title.size, 
                  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                  title.x.axis = title.x.axis,
                  title.y.axis = title.y.axis,
                  legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                  legend.text.size = legend.text.size)
  }
}


#' @title Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM for each latent process.
#' @method autoplot gmwm2
#' @param object A \code{GMWM} object.
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param color.line A \code{vector} of \code{string} that indicates the color of the lines for Empirical WV and CI. Default value is \code{c("#003C7D", "#003C7D")}
#' @param color.CI A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param line.type A \code{vector} that indicates the type of lines for Empirical WV and its CI. Default value is \code{c("solid","dotted")}
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM for each latent process.
#' @author JJB
autoplot.gmwm2 = function(object, CI = T, transparence = 0.1, color.line = c("#003C7D", "#003C7D"), 
                          color.CI = "#003C7D", line.type = c("solid","dotted"), 
                          graph.title = NA, graph.title.size= 15, 
                          axis.label.size = 13, axis.tick.size = 11, 
                          title.x.axis = expression(paste("Scale ", tau)),
                          title.y.axis = expression(paste("Wavelet Variance ", nu)),
                          legend.title = '', legend.key.size = 1, legend.title.size = 13, 
                          legend.text.size = 13, ...){
  #require package: grid
  .x=low=high=trans_breaks=trans_format=math_format=NULL
  
  # Find number of latent processes
  L = length(object$model$desc) + 1
  
  # Get names of latent processes
  nom = letters[1:L]
  
  # Construct data.frame
  df = data.frame(scales = rep(object$scales,L), WV = c(as.vector(object$decomp.theo), object$theo), 
                  process = rep(nom, each = length(object$scales)))
  WV = data.frame(var = object$wv.empir, low = object$ci.low, high = object$ci.high, scale = object$scales)
  
  p = ggplot(data = WV, mapping = aes(x = scale, y = var), colour = color.line[1]) + geom_line(linetype = line.type[1]) + geom_point(size = 3) 
    
  if(CI){
    p = p + 
      geom_line(mapping = aes(y = low), colour = color.line[2], linetype = line.type[2]) +
      geom_line(mapping = aes(y = high), colour = color.line[2], linetype = line.type[2]) +
      geom_ribbon(mapping = aes(ymin = low, ymax = high), fill = alpha(color.CI, transparence))
      #geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1, fill = "#003C7D") +
  }
  
  p = p + geom_line(data = df, mapping = aes(x = scales, y = WV, color = process)) + 
    xlab(title.x.axis ) + ylab( title.y.axis) +
    scale_y_log10( breaks = trans_breaks("log10", function(x) 10^x),
                   labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    #scale_colour_discrete(name  =" ", labels=c(object$model$desc,paste(object$model$desc, collapse = ' + '))) +
    scale_colour_discrete(name = legend.title, labels=c(object$model$desc,paste(object$model$desc, collapse = ' + '))) +
    coord_cartesian(ylim = c(min(object$ci.low), (1.05)*max(object$ci.high))) +
    xlab(title.x.axis) + ylab(title.y.axis) + ggtitle(graph.title) +
    
    theme(
      plot.title = element_text(size=graph.title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size))  
    #scale_colour_hue(name = legend.title)
  
  if (is.na(graph.title)){
    if(object$robust){
      p = p + ggtitle("Haar Wavelet Variance Robust Representation")
    }
    else{
      p = p + ggtitle("Haar Wavelet Variance Classical Representation")
    }
  }
  
  p
}


#' @title Graph Solution of the Generalized Method of Wavelet Moments
#' @description Creates a graph containing the empirical and theoretical wavelet variances constructed via GMWM.
#' @method autoplot gmwm
#' @param object A \code{GMWM} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param color.CI A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param line.type A \code{vector} that indicates the type of lines for Empirical WV, CI and Implied WV respectively. Default value is \code{c("solid","dotted","solid")}
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param facet.title.size An \code{integer} that indicates the size of facet label
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of the empirical and theoretical wavelet variance under the constructed GMWM.
#' @author JJB
#' @examples
#' # AR
#' set.seed(1336)
#' n = 200
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' mod = gmwm(2*AR1(), data = x)
#' autoplot(mod)
autoplot.gmwm = function(object,  CI = T, transparence = 0.1,  
                         color.CI = "#003C7D", line.type = c("solid","dotted","solid"), 
                         graph.title = NA, graph.title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         title.x.axis = expression(paste("Scale ", tau)),
                         title.y.axis = expression(paste("Wavelet Variance ", nu)),
                         legend.title = '', legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  #require pakage: scales, grid
  low=high=emp=theo=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
                  scale = object$scales,
                  theo = object$theo)

  p = ggplot(data = WV, mapping = aes(x = scale)) +  geom_line(aes(y = emp, colour = "LINE1"), linetype = line.type[1]) + geom_point(aes(y = emp,colour = "LINE1"), size = 3) +
    geom_line(aes(y = theo, colour = "LINE3"), linetype = line.type[3]) + 
    geom_point(aes(y = theo, colour = "LINE3"), size = 3, shape = 1) 
    
  
  if(CI){
    p = p + 
      geom_line(aes(y = low, colour = "LINE2"), linetype = line.type[2]) +
      geom_line(aes(y = high, colour = "LINE2"), linetype = line.type[2]) +
      geom_ribbon(mapping = aes(ymin = low, ymax = high), fill = alpha(color.CI, transparence))
      #geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1) 
  }
  
  #decide where to place the legend
  legendPlace = placeLegend(WV$emp[1], WV$low[ length(WV$low) ], WV$high[ length(WV$high)])
  p = p +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) 
  if (CI){
    p = p + scale_colour_manual(name=legend.title, labels=c(expression(paste("Empirical WV ", hat(nu))), 
                                                            expression(paste("CI(", hat(nu)," , 0.95)" )), expression(paste("Implied WV ", nu,"(",hat(theta),")"))), 
                                values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
      guides(colour = guide_legend(override.aes = list(size = c(1,1,1), colour = c("#000000","#999999","#56B4E9"))))
  }
  else{
    p = p + scale_colour_manual(name=legend.title, labels=c(expression(paste("Empirical WV ", hat(nu))), expression(paste("Implied WV ", nu,"(",hat(theta),")"))), 
                                values=cols, guide = guide_legend(fill = NULL,colour = NULL)) +
      guides(colour = guide_legend(override.aes = list(size = c(1,1), colour = c("#000000","#56B4E9"))))
  }
   
  p = p +
    xlab(title.x.axis) + ylab(title.y.axis) + ggtitle(graph.title) +
    theme(
          plot.title = element_text(size=graph.title.size),
          axis.title.y = element_text(size= axis.label.size),
          axis.text.y  = element_text(size= axis.tick.size),
          axis.title.x = element_text(size= axis.label.size),
          axis.text.x  = element_text(size= axis.tick.size),
          legend.key = element_rect(fill=NA), 
          legend.key.size = unit(legend.key.size, "cm"),
          legend.text = element_text(size = legend.text.size),  
          legend.title = element_text(size = legend.title.size),
          legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=legendPlace[1:2], legend.position=legendPlace[3:4])  
  #scale_colour_hue(name = legend.title)
  
  if (is.na(graph.title)){
    if(object$robust){
      p = p + ggtitle("Haar Wavelet Variance Robust Representation")
    }
    else{
      p = p + ggtitle("Haar Wavelet Variance Classical Representation")
    }
  }
  
  p
}

#' @title Compare GMWM Model Fits on Same Graph
#' @description Creates a single graph that contains two GMWM models plotted against each other.
#' @method autoplot comp
#' @param object A \code{data.frame} containing both sets of GMWM object data.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing one graphs with two GMWM models plotted against each other.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(1335)
#' n = 200
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
autoplot.comp = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#EF8A1C","LINE4"="#6E8BBF")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
                  scale = object$scales,
                  theo = object$theo1,
                  model2 = object$theo2)
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


#' @title Compare GMWM Model Fit on Split Graphs
#' @description Creates GMWM model fits of two models on split graphs within the same panel.
#' @method autoplot compSplit
#' @param object A \code{data.frame} containing both sets of GMWM object data.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing two graphs of the wavelet variance.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(1355)
#' n = 200
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = TRUE)}
autoplot.compSplit = function(object, ...){
  
  low=high=emp=theo=model2=trans_breaks=trans_format=math_format=.x=NULL
  
  cols = c("LINE1"="#000000", "LINE2"="#999999", "LINE3"="#56B4E9")
  
  WV = data.frame(emp = object$wv.empir,
                  low = object$ci.low,
                  high = object$ci.high,
                  scale = object$scales,
                  theo = object$theo1,
                  model2 = object$theo2)
  
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

#' @title Graphically Compare GMWM Model Fit Between Two Models
#' @description Creates GMWM model fits of two models on split graphs within the same panel.
#' @usage compare.models(GMWM1, GMWM2, split = FALSE)
#' @param GMWM1 A \code{gmwm} object
#' @param GMWM2 A \code{gmwm} object
#' @param split A \code{boolean} indicating true or false to place model fits on different or the same graphs.
#' @return A ggplot2 panel containing two graphs of the wavelet variance.
#' @author JJB
#' @examples
#' \dontrun{# AR
#' set.seed(8836)
#' n = 200
#' x = gen.ts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
#' GMWM1 = gmwm(AR1(), data = x)
#' GMWM2 = gmwm(2*AR1(), data = x)
#' compare.models(GMWM1, GMWM2, split = FALSE)}
compare.models = function(GMWM1, GMWM2, split = FALSE){
  x = data.frame(wv.empir = GMWM1$wv.empir, ci.low = GMWM1$ci.low, 
                 ci.high = GMWM1$ci.high, scales = GMWM1$scales, theo1 = GMWM1$theo, theo2 = GMWM2$theo) 
  if (split == TRUE){
    class(x) = "compSplit"
  }else{
    class(x) = "comp"
  }
  autoplot(x)
}
