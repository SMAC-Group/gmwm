% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Mod_cpp}
\alias{Mod_cpp}
\title{Absolute Value or Modulus of a Complex Number.}
\usage{
Mod_cpp(x)
}
\arguments{
\item{x}{A \code{cx_vec}.}
}
\value{
A \code{vec} containing the modulus for each element.
}
\description{
Computes the value of the Modulus.
}
\details{
Consider a complex number defined as: \eqn{z = x + i y} with real \eqn{x} and \eqn{y},
The modulus is defined as: \eqn{r = Mod(z) = \sqrt{(x^2 + y^2)}}
}
\examples{
Mod_cpp(c(1+.5i, 2+1i, 5+9i))
}
\keyword{internal}
