# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Maximum Overlap Discrete Wavelet Transform
#' @description Calculation of the coefficients for the discrete wavelet transformation
#' @param x A \code{vector} with dimensions N x 1. 
#' @return y A \code{field<vec>} that contains the wavelet coefficients for each decomposition level
#' @details
#' Performs a level \eqn{J} decomposition of the time series using the pyramid algorithm.
#' \eqn{J} is determined by \eqn{floor\left(log_2 \left(length\left(x\right)\right)\right)}{floor(log2(length(x)))}
#' This function was designed to minimize the amount of work a user performs. 
#' If you need more complex computations, see \code{\link{modwt_cpp}}.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' modwt(x)
modwt = function(x) {
  nlevels =  floor(log2(length(x)))
  out = .Call('gmwm_modwt_cpp', PACKAGE = 'gmwm', x, filter_name = "haar", nlevels, boundary="periodic")
  out = structure(list(data=out, nlevels=nlevels), class = "gmwm_modwt")
  invisible(out)
}

#' @title Print Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @method print gmwm_modwt
#' @export
#' @param x A \code{gmwm_modwt} object
#' @param ... further arguments passed to or from other methods.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @keywords internal
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' print(modwt(x))
print.gmwm_modwt=function(x, ...){
  cat("Results of the MODWT containing ",x$nlevels,"\n")
  print(matrix(unlist(x$data),ncol=x$nlevels,byrow=F))
}

#' @title Summary Maximum Overlap Discrete Wavelet Transform
#' @description Unlists MODWT object and places it in matrix form
#' @method summary gmwm_modwt
#' @export
#' @keywords internal
#' @param object A \code{gmwm_modwt} object
#' @param ... additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' summary(modwt(x))
summary.gmwm_modwt=function(object, ...){
  print.gmwm_modwt(object)
}