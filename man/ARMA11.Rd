% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts.model.R
\name{ARMA11}
\alias{ARMA11}
\title{Create an Autoregressive Order 1 - Moving Average Order 1 (ARMA(1,1)) Process}
\usage{
ARMA11(phi = NULL, theta = NULL, sigma2 = 1)
}
\arguments{
\item{phi}{A \code{double} containing the coefficients for \eqn{\phi _1}{phi[1]}'s for the Autoregressive 1 (AR1) term.}

\item{theta}{A \code{double} containing  the coefficients for \eqn{\theta _1}{theta[1]}'s for the Moving Average 1 (MA1) term.}

\item{sigma2}{A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of the ARMA process.}
}
\value{
An S3 object with called ts.model with the following structure:
\describe{
 \item{process.desc}{\eqn{AR1}, \eqn{MA1}, \eqn{SIGMA2}}
 \item{theta}{\eqn{\phi}{phi}, \eqn{\theta}{theta}, \eqn{\sigma^2}{sigma^2}}
 \item{plength}{Number of Parameters: 3}
 \item{obj.desc}{Depth of Parameters e.g. list(c(1,1,1))}
 \item{starting}{Guess Starting values? \code{TRUE} or \code{FALSE} (e.g. specified value)}
}
}
\description{
Sets up the necessary backend for the ARMA(1,1) process.
}
\details{
A variance is required since the model generation statements utilize 
randomization functions expecting a variance instead of a standard deviation like R.
}
\examples{
# Creates an ARMA(1,1) process with predefined coefficients.
ARMA11(phi = .23, theta = .1, sigma2 = 1)

# Creates an ARMA(1,1) process with values to be guessed on callibration.
ARMA11()
}
\author{
JJB
}
