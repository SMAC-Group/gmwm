## Creates RcppArmadillo S3 object export within R
avar = function(x, ...) UseMethod("avar")

avar.default = function(x, ...) {
  x = as.vector(x)
  av = .Call('GMWM_avar_fixed_arma', PACKAGE = 'GMWM', x)
  av$adev = sqrt(av$allan)
  av$lci = av$adev - av$errors*av$adev
  av$uci = av$adev + av$errors*av$adev
  class(av) = "avar"
  av
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
  plot(x$clusters, x$adev,log="xy",
       xlab=expression("Scale " ~ tau),
       ylab=expression("Allan Deviation " ~ phi[tau]),
       main=expression(log(tau) ~ " vs. " ~ log(phi[tau]))
  )
  lines(x$clusters, x$adev, type="l")
  lines(x$clusters, x$lci, type="l", col="grey", lty=2)
  lines(x$clusters, x$uci, type="l", col="grey", lty=2)
}

summary.avar = function(x, ...) {
  out_matrix = matrix(0, nrow = length(x$clusters), ncol = 6)
  colnames(out_matrix) = c("Time", "AVAR", "ADEV", "Lower CI", "Upper CI", "Error")
  out_matrix[,"Time"] = x$clusters
  out_matrix[,"AVAR"] = x$allan
  out_matrix[,"ADEV"] = x$adev
  out_matrix[,"Lower CI"] = x$lci
  out_matrix[,"Upper CI"] = x$uci
  out_matrix[,"Error"] = x$errors
  
  class(x) = "summary.avar"
  out_matrix
}

print.summary.avar = function(x, ...) {
  print(x)
  invisible(x)
}