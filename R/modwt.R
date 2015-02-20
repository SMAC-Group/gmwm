#' @title Maximum Overlap Discrete Wavelet Transform
#' @description Calculation of the coefficients for the discrete wavelet transformation
#' @usage modwt(x)
#' @param x A \code{vector} with dimensions N x 1. 
#' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' \eqn{J} is determined by \eqn{floor\left(log_2 \left(length(x\right)\right)\right)}{floor(log2(length(x)))}
#' This function was designed to minimize the amount of work a user performs. 
#' If you need more complex computations, see \code{\link{modwt_cpp}}.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' modwt(x)
modwt = function(x) {
  nlevels =  floor(log2(length(x)))
  out = .Call('GMWM_modwt_cpp', PACKAGE = 'GMWM', x, "haar", nlevels, "periodic")
  out = structure(list(data=out, nlevels=nlevels), class = "modwt")
  invisible(out)
}

#' @title Print Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @usage print.modwt(x, ...)
#' @param x A \code{modwt} object
#' @param ... Not used
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' print(modwt(x))
print.modwt=function(object, ...){
  cat("Results of the MODWT containing ",object$nlevels,"\n")
  print(matrix(unlist(object$data),ncol=object$nlevels,byrow=F))
}

#' @title Summary Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @usage summary.modwt(x, ...)
#' @param x A \code{modwt} object
#' @param ... Not used
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' summary(modwt(x))
summary.modwt=function(object, ...){
  print.modwt(object)
}