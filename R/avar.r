## Creates RcppArmadillo S3 object export.
avar = function(x, ...) UseMethod("avar")

avar.default <- function(x, ...) {
  x = as.vector(x)
  av = .Call('GMWM_allan_variance', PACKAGE = 'GMWM', x)
  av_out = list("clusters" = av[,1], "allan" = av[,2], "errors" = av[,3])
  class(av_out) = "avar"
  av_out
}

print.avar = function(x, ...) {
  cat("\n Clusters: \n")
  print(x$clusters, digits=5)
  cat("\n Allan Variances: \n")
  print(x$allan, digits=5)
  cat("\n Errors: \n")
  print(x$errors, digits=5)
}

plot.avar = function(x, ...){
  
  
}

summary.avar = function(x, ...) {
  class(x) = "summary.avar"
  x
}

print.summary.avar = function(x, ...) {
  print.avar(x)
  invisible(x)
}