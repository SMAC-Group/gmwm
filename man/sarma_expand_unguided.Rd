% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sarma_expand_unguided}
\alias{sarma_expand_unguided}
\title{(Internal) Expand the SARMA Parameters}
\usage{
sarma_expand_unguided(params, np, nq, nsp, nsq, ns, p, q)
}
\arguments{
\item{params}{A \code{vec} containing the theta values of the parameters.}

\item{np}{An \code{unsigned int} containing the number of non-seasonal phi parameters.}

\item{nq}{An \code{unsigned int} containing the number of non-seasonal theta parameters.}

\item{nsp}{An \code{unsigned int} containing the number of seasonal phi parameters.}

\item{nsq}{An \code{unsigned int} containing the number of seasonal theta parameters.}

\item{p}{An \code{unsigned int} that is the total size of the phi vector.}

\item{q}{An \code{unsigned int} that is the total size of the theta vector.}
}
\value{
A \code{field<vec>} that contains the expansion.
}
\description{
(Internal) Expand the SARMA Parameters
}
\keyword{internal}
