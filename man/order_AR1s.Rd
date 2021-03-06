% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{order_AR1s}
\alias{order_AR1s}
\title{Order AR1s by size of phi.}
\usage{
order_AR1s(theta, desc, objdesc)
}
\arguments{
\item{theta}{A \code{vec} with dimensions N x 1 that contains user-supplied initial values for parameters}

\item{desc}{A \code{vector<string>} indicating the models that should be considered.}

\item{objdesc}{A \code{field<vec>} containing a list of parameters (e.g. AR(1) = c(1,1), ARMA(p,q) = c(p,q,1))}
}
\value{
A \code{vec} that has AR1s shown in descending parameter value.
}
\description{
Changes the order of AR1s in a string by size.
}
\keyword{internal}
