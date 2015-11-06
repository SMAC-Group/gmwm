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

#' @title Calculate the Allan Variance
#' @description Computes the Allan Variance
#' @param x A \code{vec} containing the time series under observation.
#' @return Allan variance fixed
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' avar(x)
avar = function(x) {
  x = as.vector(x)
  av = .Call('gmwm_avar_mo_cpp', PACKAGE = 'gmwm', x)
  av = list(clusters = av[,1], allan=av[,2], errors=av[,3])
  av$adev = sqrt(av$allan)
  av$lci = av$adev - av$errors*av$adev
  av$uci = av$adev + av$errors*av$adev
  class(av) = "avar"
  av
}

#' @title Prints Allan Variance
#' @description Displays the allan variance information
#' @method print avar
#' @export
#' @param x A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @author JJB
#' @return console output
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' print( out )
print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

#' @title Plot Allan Variance
#' @description Displays a plot containing the allan variance
#' @method plot avar
#' @export
#' @param x A \code{avar} object.
#' @param ... Arguments to be passed to methods
#' @author JJB
#' @return ggplot2 graph
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' plot( out )
plot.avar = function(x, ...){
  plot(x$clusters, x$adev,log="xy",
       xlab=expression("Scale " ~ tau),
       ylab=expression("Allan Deviation " ~ phi[tau]),
       main=expression(log(tau) ~ " vs. " ~ log(phi[tau]))
  )
  lines(x$clusters, x$adev, type="l")
  lines(x$clusters, x$lci, type="l", col="grey", lty=2)
  lines(x$clusters, x$uci, type="l", col="grey", lty=2)
}


#' @title Summary Allan Variance
#' @description Displays the summary table of allan variance
#' @method summary avar
#' @export
#' @param object A \code{avar} object.
#' @param ...  additional arguments affecting the summary produced.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = avar(x)
#' summary( out )
summary.avar = function(object, ...) {
  out_matrix = matrix(0, nrow = length(object$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "AVAR", "ADEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = object$clusters
  out_matrix[,"AVAR"] = object$allan
  out_matrix[,"ADEV"] = object$adev
  out_matrix[,"Lower CI"] = object$lci
  out_matrix[,"Upper CI"] = object$uci
  out_matrix[,"Error"] = object$errors
  
  class(out_matrix) = "summary.avar"
  out_matrix
}