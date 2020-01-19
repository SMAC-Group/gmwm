#' @title Calculate the Hadamard Variance
#' @description Computes the Hadamard Variance
#' @param x     A \code{vec} containing the time series under observation.
#' @param type  A \code{string} containing either \code{"mo"} for Maximal Overlap or \code{"to"} for Tau Overlap
#' @return Hadamard variance fixed
#' @return hadam A \code{list} that contains:
#' \itemize{
#'  \item{"clusters"}{The size of the cluster}
#'  \item{"hadamard"}{The Hadamard variance}
#'  \item{"errors"}{The error associated with the variance estimation.}
#' }
#' @details 
#' The decomposition and the amount of time it takes to perform it depends on whether you are using
#' the Tau Overlap or the Maximal Overlap.
#' 
#' Maximal Overlap Hadamard Variance
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/3}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log3(N))}}
#' Then, \eqn{M = N - 3n} samples exist. 
#' The Maximal-overlap estimator is given by:
#' \eqn{\frac{1}{{6(M - 3m + 1)}}\sum\limits_{i = 1}^{M - 3m + 1} {{{[{y_{i + 2}} - 2{y_{i + 1}} + {y_i}]}^2}}}{See PDF Manual}
#' 
#' where \eqn{{y_i} = {{\bar y}_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}} .}{ See PDF Manual}.
#' 
#' Tau-Overlap Hadamard Variance
#' Given \eqn{N} equally spaced samples with averaging time \eqn{\tau = n\tau _0}{tau = n*tau_0},
#' where \eqn{n} is an integer such that \eqn{ 1 \le n \le \frac{N}{2}}{1<= n <= N/2}.
#' Therefore, \eqn{n} is able to be selected from \eqn{\left\{ {n|n < \left\lfloor {{{\log }_2}\left( N \right)} \right\rfloor } \right\}}{{n|n< floor(log2(N))}}
#' Then, a sampling of \eqn{m = \left\lfloor {\frac{{N - 1}}{n}} \right\rfloor  - 1} samples exist. 
#' The tau-overlap estimator is given by:
#' \eqn{\frac{1}{{6(M - 2)}}\sum\limits_{t = 1}^{M - 2} {{{\left( {{y_{t + 2}} - 2{y_{t + 1}} + {y_t}} \right)}^2}}}{See PDF Manual}
#' where \eqn{ {\bar y_t}\left( \tau  \right) = \frac{1}{\tau }\sum\limits_{i = 0}^{\tau  - 1} {{{\bar y}_{t - i}}}  }{See PDF Manual}.
#' 
#' 
#' @author Avinash Balakrishnan, JJB
#' @examples
#' set.seed(999)
#' # Simulate white noise (P 1) with sigma^2 = 4
#' N = 100000
#' white.noise = rnorm(N, 0, 2)
#' #plot(white.noise,ylab="Simulated white noise process",xlab="Time",type="o")
#' #Simulate random walk (P 4)
#' random.walk = cumsum(0.1*rnorm(N, 0, 2))
#' combined.ts = white.noise+random.walk
#' hadam_mo = hadam(combined.ts)
#' 
#' hadam_to = hadam(combined.ts, type = "to")
hadam = function(x, type = "mo") {
  x = as.vector(x)
  
  if(type == "mo"){
    hv = .Call('_gmwm_hadam_mo_cpp', PACKAGE = 'gmwm', x)
  }else{
    hv = .Call('_gmwm_hadam_to_cpp', PACKAGE = 'gmwm', x)
  }
  
  hv = list(clusters = hv[,1], hadamard=hv[,2], errors=hv[,3])
  hv$hdev = sqrt(hv$hadamard)
  hv$lci = hv$hdev - hv$errors*hv$hdev
  hv$uci = hv$hdev + hv$errors*hv$hdev
  hv$type = type
  class(hv) = "hadam"
  hv
}

#' @title Prints Hadamard Variance
#' @description Displays the all Hadamard variance information
#' @method print hadam
#' @export
#' @param x   A \code{hadam} object.
#' @param ... Arguments to be passed to methods
#' @author JJB
#' @return console output
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = hadam(x)
#' print( out )
print.hadam = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Hadamard Variances: \n")
  print(x$hadamard, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' @title Wrapper to ggplot Hadamard Variance Graph
#' @description Displays a plot containing the Hadamard variance
#' @method plot hadam
#' @export
#' @param x   A \code{hadam} object.
#' @template CommonParams
#' @author JJB, Wenchao
#' @return A ggplot2 graph containing the Hadamard variance.
#' @note 
#' It is known that the calculation for confidence interval is incorrect, therefore \code{CI}
#' can only be set to FALSE currently.
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = hadam(x)
#' plot( out )
plot.hadam = function(x, CI = F, transparence = 0.1, background = 'white', bw = F, 
                      CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                      point.size = NULL, point.shape = NULL,
                      title = NULL, title.size= 15, 
                      axis.label.size = 13, axis.tick.size = 11, 
                      axis.x.label = expression("Cluster "~tau~"(sec)"),
                      axis.y.label = expression("Hadamard Variance "~phi[tau]),
                      legend.title = '',  legend.label = NULL,
                      legend.key.size = 1, legend.title.size = 13, 
                      legend.text.size = 13, ...){
  
  autoplot.hadam(x, CI = CI, transparence = transparence, background = background, bw = bw, 
                 CI.color = CI.color, line.type = line.type, line.color = line.color,
                 point.size = point.size, point.shape = point.shape,
                 title = title, title.size= title.size, 
                 axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                 axis.x.label = axis.x.label,
                 axis.y.label = axis.y.label,
                 legend.title = legend.title, legend.label = legend.label,
                 legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                 legend.text.size = legend.text.size )
  
}

#' @title Graph Hadamard Variance
#' @description Displays a plot containing the Hadamard variance
#' @method autoplot hadam
#' @export
#' @param object A \code{hadam} object.
#' @template CommonParams
#' @author JJB, Wenchao
#' @return A ggplot2 graph containing the Hadamard variance.
#' @note 
#' It is known that the calculation for confidence interval is incorrect, therefore \code{CI}
#' can only be set to FALSE currently.
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = hadam(x)
#' autoplot( out )
autoplot.hadam = function(object, CI = F, transparence = 0.1, background = 'white', bw = F, 
                          CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                          point.size = NULL, point.shape = NULL,
                          title = NULL, title.size= 15, 
                          axis.label.size = 13, axis.tick.size = 11, 
                          axis.x.label = expression("Cluster "~tau~"(sec)"),
                          axis.y.label = expression("Hadamard Variance "~phi[tau]),
                          legend.title = '',  legend.label =  NULL,
                          legend.key.size = 1, legend.title.size = 13, 
                          legend.text.size = 13, ...){
  if(CI == T){
    warning("It is known that the calculation for confidence interval is incorrect, therefore 'CI' can only be set to FALSE currently.")
    CI = F
  }
  
  # The format of object that can be passed to graphingVar()
  names(object) = c("scales", "variance", "errors", "hdev", "ci_low", "ci_high", "type")
  
  #check parameter
  params = 'legend.label'
  if(CI){
    #Update later
    #requireLength = 2
    #legend.label.default = c(bquote("Empirical AV"~hat(sigma)^2), bquote("CI("*hat(nu)*", "*.(1 - object$alpha)*")" ))
  }else{
    requireLength = 1
    legend.label.default = c(bquote("Empirical HV "~hat(sigma)^2))
  }
  
  default = list(legend.label.default)
  nullIsFine = T
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  p = graphingVar(object, CI = CI, transparence = transparence, background = background, bw = bw, 
                  CI.color = CI.color, line.type = line.type, line.color = line.color,
                  point.size = point.size, point.shape = point.shape,
                  title = title, title.size= title.size, 
                  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                  axis.x.label = axis.x.label,
                  axis.y.label = axis.y.label,
                  legend.title = legend.title,  legend.label = legend.label,
                  legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                  legend.text.size = legend.text.size )
  
  
  if(is.null(title)){
    
    if (object$type == "mo"){
      type = "Maximal Overlap"
    }else{
      type = "Tau Overlap"
    }
    
    p = p +
      ggtitle(paste0("Hadamard Variance (",type, ")"))
    
  }
  
  p
  
}

#' @title Summary Hadamard Variance
#' @description Displays the summary table of Hadamard variance
#' @method summary hadam
#' @export
#' @param object A \code{hadam} object.
#' @param ...    Additional arguments affecting the summary produced.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = hadam(x)
#' summary( out )
summary.hadam = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "hadam", "HDEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = object$clusters
  out_matrix[,"hadam"] = object$hadamard
  out_matrix[,"HDEV"] = object$hdev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci
  out_matrix[,"Error"] = object$errors
  
  class(out_matrix) = "summary.hadam"
  out_matrix
}