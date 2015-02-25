#' @title Calculate the Asymptotic Covariance Matrix
#' @description Places the Asymptotic Covariance Matrix in print form.
#' @usage wvcov(signal.modwt, signal.wvar, compute.v="diag")
#' @param signal.modwt A \code{modwt} object that contains the modwt decomposition.
#' @param signal.wvar A \code{wvar} object that contains the wavelet variance.
#' @param compute.v A \code{string} that indicates the type of covariance matrix to compute. Supports: "diag"
#' @return A \code{list} with the structure:
#' \itemize{
#'   \item{"V"}{Covariance Matrix}
#'   \item{"V.r"}{Covariance Matrix}
#'   \item{"nlevels"}{Level of decomposition J}
#'   \item{"compute.v"}{Type of Covariance Matrix}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"scales"}{Tau scales (2^(1:J))}
#'   \item{"wv.empir"}{Empirical Wavelet Variance}
#' }
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(decomp)
#' out = wvcov(decomp, wv, compute.v="diag")
#' 
#' # Robust
#' decomp = modwt(x)
#' wv = wvar(decomp, robust = TRUE)
#' out = wvcov(decomp, wv, compute.v="diag")
wvcov = function(signal.modwt, signal.wvar, compute.v="diag"){
  
  if(!is(signal.modwt,"modwt")){
    stop("Need to supply a modwt object as the first parameter.")
  }
  
  if(!is(signal.wvar,"wvar")){
    stop("Need to supply a wvar object as the second parameter.")
  }
  
  out = .Call('GMWM_compute_cov_cpp', PACKAGE = 'GMWM', signal.modwt$data, signal.modwt$nlevels, compute.v, signal.wvar$robust, signal.wvar$eff)
  out = structure(list(V=out[[1]],
                       V.robust=out[[2]], 
                       nlevels=signal.modwt$nlevels, 
                       compute.v = compute.v, 
                       robust = signal.wvar$robust, 
                       eff = signal.wvar$eff, 
                       scales = signal.wvar$scales,
                       wv.empir = signal.wvar$variance,
                       ci_low = signal.wvar$ci_low,
                       ci_high = signal.wvar$ci_high), class = "wvcov")
  invisible(out)
}



#' @title Print Asymptotic Covariance Matrix
#' @description Places the Asymptotic Covariance Matrix in print form.
#' @method print wvcov
#' @param x A \code{wvcov} object
#' @param ...  further arguments passed to or from other methods
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(decomp)
#' print(wvcov(decomp,wv,compute.v="diag"))
print.wvcov = function(x, ...){
  print(x$V)
}

#' @title Summary Wavelet Covariance Matrix
#' @description Prints the Wavelet Covariance Matrix
#' @method summary wvcov
#' @param object A \code{wvcov} object
#' @param ...  additional arguments affecting the summary produced.
#' @return Prints the modwt matrix decomposition
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' decomp = modwt(x)
#' wv = wvar(decomp)
#' summary(wvcov(decomp,wv,compute.v="diag"))
summary.wvcov=function(object, ...){
  
  name = if(object$robust){
    "robust" 
  }else{
    "classical"
  }
  cat("The asymptotic ", object$compute.v, " using the ",name, " method.\n",sep="")
  if(object$robust){
    cat("The robust form was created using efficiency=",object$eff,"\n",sep="")
  }  
  
  print(object)
}