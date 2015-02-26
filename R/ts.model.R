#' @title Integer Check
#' @description Checks whether the submitted value is an integer
#' @usage is.whole(x)
#' @param x A \code{numeric} value to check to see if it is an integer.
#' @return A \code{boolean} value indicating whether the value is an integer or not.
#' @author JJB
#' @examples
#' is.whole(2.3)
#' is.whole(4)
#' is.whole(c(1,2,3))
#' is.whole(c(.4,.5,.6))
#' is.whole(c(7,.8,9))
is.whole = function(x){ is.numeric(x) && floor(x)==x } 

#' @title Create an Autoregressive 1 [AR(1)] Process
#' @description Setups the necessary backend for the AR1 process.
#' @usage AR1(phi = 0.91, sigma = 1)
#' @param phi A \code{double} value for the \eqn{\phi}{phi} of an AR1 process.
#' @param sigma A \code{double} value for the \eqn{\sigma}{sigma} process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"AR1"}
#'  \item{theta}{\eqn{phi}{\phi}, \eqn{sigma}{\sigma}}
#' }
#' @author JJB
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.

#' @examples
#' AR1()
#' AR1(phi=.32, sigma=1.3)
AR1 = function(phi = 0.91, sigma = 1) {
  if(length(phi) != 1 & length(sigma) != 1){
    stop("Bad AR1 model submitted. Must be double values for two parameters.")
  }
  out = structure(list(desc = "AR1",
                       theta = c(phi,sigma),
                       plength = 2), class = "ts.model")
  invisible(out)
}

#' @title Create an Quantisation Noise (QN) Process
#' @description Sets up the necessary backend for the QN process.
#' @usage QN(q2)
#' @param q2 A \code{double} value for the \eqn{Q^2}{Q^2} of a QN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"QN"}
#'  \item{theta}{\eqn{Q^2}{Q^2}}
#' }
#' @author JJB
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.

#' @examples
#' QN()
#' QN(q2=3.4)
QN = function(q2 = 1) {
  if(length(q2) != 1){
    stop("Bad QN model submitted. Must be a double that indicates the Q2 value.")
  }
  out = structure(list(desc = "QN",
                       theta = q2,
                       plength = 1), class = "ts.model")
  invisible(out)
}

#' @title Create an White Noise (WN) Process
#' @description Sets up the necessary backend for the WN process.
#' @usage WN(sigma=1)
#' @param sigma A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"WN"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#' }
#' @author JJB
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.

#' @examples
#' WN()
#' WN(sigma=3.4)
WN = function(sigma = 1) {
  if(length(sigma) != 1){
    stop("Bad WN model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(desc = "WN",
                       theta = sigma,
                       plength = 1), class = "ts.model")
  invisible(out)
}

#' @title Create an Random Walk (RW) Process
#' @description Sets up the necessary backend for the RW process.
#' @usage RW(sigma=1)
#' @param sigma A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of a RW process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"RW"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#' }
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.
#' @author JJB

#' @examples
#' RW()
#' RW(sigma=3.4)
RW = function(sigma = 1) {
  if(length(sigma) != 1){
    stop("Bad RW model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(desc = "RW",
                       theta = sigma,
                       plength = 1), class = "ts.model")
  invisible(out)
}

#' @title Create an Drift (DR) Process
#' @description Sets up the necessary backend for the DR process.
#' @usage DR(slope = 5)
#' @param slope A \code{double} value for the slope of a DR process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"DR"}
#'  \item{theta}{slope}
#' }
#' @author JJB

#' @examples
#' DR()
#' DR(slope=3.4)
DR = function(slope = 5) {
  if(length(slope) != 1){
    stop("Bad Drift model submitted. Must be a double that indicates a slope.")
  }
  out = structure(list(desc = "DR",
                       theta = slope,
                       plength = 1), class = "ts.model")
  invisible(out)
}

#' @title Create an Autoregressive Moving Average (ARMA) Process
#' @description Sets up the necessary backend for the ARMA process.
#' @usage ARMA(ar = 1, ma = 1, sigma = 1.0)
#' @param ar A \code{vector} or \code{integer} containing either the coefficients for \eqn{\phi}{phi}'s or the process number \eqn{p} for the Autoregressive (AR) term.
#' @param ma A \code{vector} or \code{integer} containing either the coefficients for \eqn{\theta}{theta}'s or the process number \eqn{q} for the Moving Average (MA) term.
#' @param sigma A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of the ARMA process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{\eqn{AR*p}{AR x p}, \eqn{MA*q}{MA x q}}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#' }
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.
#' @author JJB

#' @examples
#' # Create an ARMA(1,2) process
#' ARMA(ar=1,2)
#' # Creates an ARMA(3,2) process with predefined coefficients.
#' ARMA(ar=c(0.23,.43, .59), ma=c(0.4,.3))
#' 
#' # Creates an ARMA(3,2) process with predefined coefficients and standard deviation
#' ARMA(ar=c(0.23,.43, .59), ma=c(0.4,.3), sigma = 1.5)
ARMA = function(ar = 1, ma = 1, sigma = 1.0) {
  if(length(ar) == 1 & length(ma) == 1){
    if(is.whole(ar) & is.whole(ma)){
      ar = rep(0, ar)
      ma = rep(1, ma)
    }
  }
  out = structure(list(desc = c(rep("AR", length(ar)), rep("MA",length(ma)), "SIGMA2"),
                       theta = c(ar, ma, sigma),
                       plength = length(ar)+length(ma) + 1), class = "ts.model")
  invisible(out)
}

#' @title Multiple a ts.model by constant
#' @description Sets up the necessary backend for creating multiple model objects.
#' @method * ts.model
#' @param x A \code{numeric} value
#' @param y A \code{ts.model} object
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{y desc replicated x times}
#'  \item{theta}{y theta replicated x times}
#' }
#' @author JJB
#' @export
#' @examples
#' 4*DR()+2*WN()
#' DR()*4 + WN()*2
#' AR1(phi=.3,sigma=.2)*3
"*.ts.model" = function(x, y) {
  # Handles the ts.model()*c case
  if(!is.numeric(x)){
    temp = x
    x = y
    y = temp
  }
  out = structure(list(desc = rep(y$desc,x),
                       theta = rep(y$theta,x),
                       plength = rep(y$plength,x)), class = "ts.model")
  invisible(out)
}

#' @title Add ts.model objects together
#' @description Sets up the necessary backend for combining ts.model objects.
#' @method + ts.model
#' @param x A \code{ts.model} object
#' @param y A \code{ts.model} object
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{combined x, y desc}
#'  \item{theta}{combined x, y theta}
#' }
#' @author JJB
#' @export
#' @examples
#' DR()+WN()
#' AR1(phi=.3,sigma=.2)
"+.ts.model" = function(x, y) {
  out = structure(list(desc = c(x$desc, y$desc),
                       theta = c(x$theta,y$theta),
                       plength = c(x$plength, y$plength)), class = "ts.model")
  invisible(out)
}

#' @title Multiple a ts.model by constant
#' @description Sets up the necessary backend for creating multiple model objects.
#' @method print ts.model
#' @param x A \code{numeric} value
#' @param ... further arguments passed to or from other methods.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}
#'  \item{theta}
#' }
#' @author JJB
#' @examples
#' DR()+WN()+RW()+AR1()+ARMA(1,2)
print.ts.model = function(x, ...) cat("Desc:\n", x$desc, "\nValues:\n", x$theta)
