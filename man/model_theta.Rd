% Generated by roxygen2: do not edit by hand
% Please edit documentation in src/ts_model_cpp.cpp, src/ts_model_cpp.h
\name{model_theta}
\alias{model_theta}
\title{Generate the ts model object's theta vector}
\usage{
model_theta(desc)
}
\arguments{
\item{desc}{A \code{vector<string>} that contains a list of the strings of each process.}
}
\value{
A \code{vec} with values initialized at 0 that span the space of parameters to be estimated.
}
\description{
Creates the ts.model's theta vector
}
\details{
This function currently does NOT support ARMA(P,Q) models. 
That is, there is no support for ARMA(P,Q), AR(P), or MA(Q).
There is support for ARMA11, AR1, MA1, GM, WN, DR, QN, and RW.
}
\keyword{internal}
