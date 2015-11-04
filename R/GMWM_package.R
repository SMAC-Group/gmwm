#' Generalized Method of Wavelet Moments (GMWM) Package
#'
#' Generalized Method of Wavelet Moments (GMWM) is an estimation technique for the parameters of time series models. It uses the wavelet variance in a moment matching approach that makes it particularly suitable for the estimation of certain state-space models. Furthermore, there exists a robust implementation of GMWM, which allows the robust estimation of some state-space models and ARIMA models. Lastly, the package provides the ability to quickly generate time series data, perform different wavelet decompositions, and visualizations. 
#' 
#' @details 
#' \tabular{ll}{
#' Package: \tab GMWM\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0.0\cr
#' Date: \tab 2015-11-02\cr
#' License: \tab file LICENSE\cr
#' }
#' 
#' @author
#' James Balamuta \email{balamut2@@illinois.edu},
#' Stephane Guerrier \email{stephane@@illinois.edu},
#' Roberto Molinari \email{roberto.molinari@@unige.ch},
#' Wenchao Yang \email{wyang40@@illinois.edu}
#'
#' Stephane Guerrier \email{stephane@@illinois.edu}
#' @name gmwm-package
#' @docType package
#' @useDynLib gmwm
#' @importFrom Rcpp evalCpp
#' @importFrom devtools install_github
#' @importFrom grDevices gray.colors hcl
#' @importFrom graphics lines plot
#' @importFrom methods is
#' @importFrom stats arima predict
#' @importFrom utils install.packages installed.packages tail
#' @import ggplot2 grid scales reshape2 gridExtra
#' @exportPattern ^[[:alpha:]]+
NULL
