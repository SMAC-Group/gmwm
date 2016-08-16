# Copyright (C) 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
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


#' @title Auto-Covariance and Correlation Functions
#' @description The acf function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @param x      A \code{matrix} with dimensions \eqn{N \times S}{N x S} or N observations and S processes
#' @param lagmax A \code{integer} indicating the max lag.
#' @param cor    A \code{bool} indicating whether the correlation 
#' (\code{TRUE}) or covariance (\code{FALSE}) should be computed.
#' @param demean A \code{bool} indicating whether the data should be detrended
#'  (\code{TRUE}) or not (\code{FALSE})
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @details 
#' \code{lagmax} default is \eqn{10*log10(N/m)} where \eqn{N} is the number of
#' observations and \eqn{m} is the number of series being compared. If 
#' \code{lagmax} supplied is greater than the number of observations, then one
#' less than the total will be taken.
#' @examples 
#' # Get Autocorrelation
#' m = ACF(datasets::AirPassengers)
#' 
#' # Get Autocovariance and do not remove trend from signal
#' m = ACF(datasets::AirPassengers, cor = FALSE, demean = FALSE)
ACF = function(x, lagmax = 0, cor = TRUE, demean = TRUE){

  # Force to matrix form
  if(is.ts(x) || is.atomic(x)){
    x2 = data.matrix(x)        
  }
    
  o = .acf(x2, lagmax, cor, demean)
  
  # Try to get names...
  name_vars = colnames(x)
  
  # No names? No sweat. Just use NSE!
  if(is.null(name_vars)){
    name_vars = deparse(substitute(x))
  }
  
  # Make pretty names
  dimnames(o)  = list(seq_len(nrow(o)), name_vars, name_vars)
  
  o = structure(o, type = (cor == 1), n = nrow(x2), class = c("ACF","array"))
  o
}


#' @title Auto-Covariance and Correlation Functions
#' @description The acf function computes the estimated
#' autocovariance or autocorrelation for both univariate and multivariate cases.
#' @param x,object  An \code{"ACF"} object from \code{\link{ACF}}.
#' @param show.ci   A \code{bool} indicating whether to show confidence region
#' @param ci        A \code{double} containing the 1-alpha level. Default is 0.95
#' @param ...       Additional parameters
#' @return An \code{array} of dimensions \eqn{N \times S \times S}{N x S x S}.
#' @rdname plot.ACF
#' @export
#' @examples 
#' # Calculate the Autocorrelation
#' m = ACF(datasets::AirPassengers)
#' 
#' # Plot with 95% CI
#' plot(m) 
#' 
#' # Plot with 90% CI
#' plot(m, ci = 0.90) 
#' 
#' # Plot without 95% CI
#' plot(m, show.ci = FALSE)
plot.ACF = function(x, show.ci = TRUE, ci = 0.95, ...){
  autoplot.ACF(object = x, show.ci = show.ci, ci = ci, ...)
}

#' @rdname plot.ACF
#' @export
autoplot.ACF = function(object, show.ci = TRUE, ci = 0.95, ...){
  
  # Quiet the warnings...
  Lag = xmin = xmax = ymin = ymax = NULL 
  
  # Wide to long array transform
  x2 = as.data.frame.table(object, responseName = "ACF")
  
  colnames(x2) = c("Lag", "Signal X", "Signal Y", "ACF")
  
  # Remove character cast
  x2$Lag = as.numeric(x2$Lag)
  
  # Create key
  x2$key = ifelse(
    x2$`Signal X` == x2$`Signal Y`,
    paste0(x2$`Signal X`),
    paste0(x2$`Signal X`, " & ", x2$`Signal Y`)
  )
  
  # Plot with facetting
  g = ggplot(data = x2, mapping = aes(x = Lag, y = ACF)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = Lag, yend = 0)) +
    facet_wrap(~ key) + theme_bw()
  
  
  if(show.ci){
    
    clim0 = qnorm( (1 + ci)/2 ) / sqrt(attr(object,'n'))

    ci.region = data.frame(xmin=-Inf, xmax=Inf, ymin=-clim0, ymax=clim0)

    g = g + geom_rect(data = ci.data, 
                      aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
                      fill = "blue", alpha = 0.10,
                      inherit.aes = FALSE)
  }
  
  if(attr(object,'type')){
    g = g + ylab("Autocorrelation")
  } else{
    g = g + ylab("Autocovariance")
  }
  
  g
}
