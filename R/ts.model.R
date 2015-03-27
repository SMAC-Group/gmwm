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
#' @param phi A \code{double} value for the \eqn{\phi}{phi} of an AR1 process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"AR1"}
#'  \item{theta}{\eqn{\phi}{phi}, \eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1,1)}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' AR1()
#' AR1(phi=.32, sigma=1.3)
AR1 = function(phi = NULL, sigma2 = NULL) {
  adv = TRUE;
  if(is.null(phi) || is.null(sigma2)){
    phi = 0;
    sigma2 = 1;
    adv = FALSE;
  }
  if(length(phi) != 1 & length(sigma2) != 1){
    stop("Bad AR1 model submitted. Must be double values for two parameters.")
  }
  out = structure(list(desc = c("AR1","SIGMA2"),
                       theta = c(phi,sigma2),
                       plength = 2,
                       obj.desc = "AR1",
                       obj = list(c(1,1)),
                       adv = adv), class = "ts.model")
  invisible(out)
}

#' @title Create an Quantisation Noise (QN) Process
#' @description Sets up the necessary backend for the QN process.
#' @param q2 A \code{double} value for the \eqn{Q^2}{Q^2} of a QN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"QN"}
#'  \item{theta}{\eqn{Q^2}{Q^2}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1)}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' QN()
#' QN(q2=3.4)
QN = function(q2 = NULL) {
  adv=TRUE
  if(is.null(q2)){
    q2 = 2
    adv = FALSE
  }
  if(length(q2) != 1){
    stop("Bad QN model submitted. Must be a double that indicates the Q2 value.")
  }
  out = structure(list(desc = "QN",
                       theta = q2,
                       plength = 1,
                       obj.desc = "QN",
                       obj = list(1),
                       adv = adv), class = "ts.model")
  invisible(out)
}

#' @title Create an White Noise (WN) Process
#' @description Sets up the necessary backend for the WN process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"WN"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1)}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' WN()
#' WN(sigma=3.4)
WN = function(sigma2 = NULL) {
  adv=TRUE
  if(is.null(sigma2)){
    sigma2 = 3
    adv = FALSE
  }
  if(length(sigma2) != 1){
    stop("Bad WN model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(desc = "WN",
                       theta = sigma2,
                       plength = 1,
                       obj.desc = "WN",
                       obj = list(1),
                       adv = adv), class = "ts.model")
  invisible(out)
}

#' @title Create an Random Walk (RW) Process
#' @description Sets up the necessary backend for the RW process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"RW"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1)}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' RW()
#' RW(sigma=3.4)
RW = function(sigma2 = NULL) {
  adv=TRUE
  if(is.null(sigma2)){
    sigma2 = 4
    adv = FALSE
  }
  if(length(sigma2) != 1){
    stop("Bad RW model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(desc = "RW",
                       theta = sigma2,
                       plength = 1,
                       obj.desc = "RW",
                       obj = list(1),
                       adv = adv), class = "ts.model")
  invisible(out)
}

#' @title Create an Drift (DR) Process
#' @description Sets up the necessary backend for the DR process.
#' @param slope A \code{double} value for the slope of a DR process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{"DR"}
#'  \item{theta}{slope}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1)}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' DR()
#' DR(slope=3.4)
DR = function(slope = NULL) {
  adv=TRUE
  if(is.null(slope)){
    slope = 5
    adv = FALSE
  }
  if(length(slope) != 1){
    stop("Bad Drift model submitted. Must be a double that indicates a slope.")
  }
  out = structure(list(desc = "DR",
                       theta = slope,
                       plength = 1,
                       obj.desc = "DR",
                       obj = list(1),
                       adv = adv), class = "ts.model")
  invisible(out)
}

#' @title Create an Autoregressive Moving Average (ARMA) Process
#' @description Sets up the necessary backend for the ARMA process.
#' @param ar A \code{vector} or \code{integer} containing either the coefficients for \eqn{\phi}{phi}'s or the process number \eqn{p} for the Autoregressive (AR) term.
#' @param ma A \code{vector} or \code{integer} containing either the coefficients for \eqn{\theta}{theta}'s or the process number \eqn{q} for the Moving Average (MA) term.
#' @param sigma2 A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of the ARMA process.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}{\eqn{AR*p}{AR x p}, \eqn{MA*q}{MA x q}}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(c(length(ar),length(ma),1) )}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
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
#' ARMA(ar=c(0.23,.43, .59), ma=c(0.4,.3), sigma2 = 1.5)
ARMA = function(ar = 1, ma = 1, sigma2 = 1.0) {
  adv = TRUE
  if(length(ar) == 1 & length(ma) == 1){
    if(is.whole(ar) & is.whole(ma)){
      ar = rep(-1, ar)
      ma = rep(-2, ma)
      adv = FALSE
    }
  }
  out = structure(list(desc = c(rep("AR", length(ar)), rep("MA",length(ma)), "SIGMA2"),
                       theta = c(ar, ma, sigma2),
                       plength = length(ar)+length(ma) + 1,
                       obj.desc = "ARMA",
                       obj = list(c(length(ar),length(ma),1)),
                       adv = adv), class = "ts.model")
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
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(c(1,1),c(1,1))}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
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
                       plength = y$plength*x,
                       obj.desc = rep(y$obj.desc,x),
                       obj = rep(y$obj,x),
                       adv = y$adv), class = "ts.model")
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
#'  \item{plength}{Number of Parameters}
#'  \item{obj}{Depth of Parameters e.g. list(1, c(1,1), c(length(ar),length(ma),1) )}
#'  \item{adv}{Advanced Input TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @export
#' @examples
#' DR()+WN()
#' AR1(phi=.3,sigma=.2)
"+.ts.model" = function(x, y) {
  adv = FALSE
  if(y$adv & x$adv){
    adv = TRUE
  }
  out = structure(list(desc = c(x$desc, y$desc),
                       theta = c(x$theta,y$theta),
                       plength = x$plength + y$plength,
                       obj.desc = c(x$obj.desc, y$obj.desc),
                       obj = c(x$obj, y$obj),
                       adv = adv), class = "ts.model")
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
print.ts.model = function(x, ...){

  desctable = data.frame("Terms" = x$desc, "Starting Values" = x$theta);
  cat("\nAdvanced Enabled:", x$adv, "\n")
  if(!x$adv){
    cat("The program will attempt to guess starting values for...\n")
    print(desctable[,1], row.names = FALSE)
    cat("This model will only work using the gmwm() function.\n",
        "To have the option of using adv.gmwm(), please supply values for each parameter.\n")
  }else{
    print(desctable, row.names = FALSE)
    cat("This model is able to be supplied to both gmwm() or adv.gmwm()\n")
  }
}