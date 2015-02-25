#' @title Wavelet Variance
#' @description Calculates the (MODWT) wavelet variance
#' @usage wvar(x, p = 0.025, robust = FALSE, eff = 0.6)
#' @param x A \code{modwt} object that contains the modwt decomposition.
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param p A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @return A \code{list} with the structure:
#' \itemize{
#'   \item{"variance"}{Wavelet Variance},
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"p"}{p value used for CI}
#' }
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' # Default
#' wvar(modwt(x))
#' # Robust
#' wvar(modwt(x), robust = TRUE, eff=0.3)
#' # Different p/2 value (e.g. alpha = 0.1)
#' wvar(modwt(x), p = 0.05)
wvar = function(x, p = 0.025, robust = FALSE, eff = 0.6) {
  if(!is(x,"modwt")){
    stop("Need to supply the modwt class.")
  }
  out = .Call('GMWM_wvar_cpp', PACKAGE = 'GMWM', x$data, x$nlevels, robust, eff, p, "eta3", "haar")
  scales = .Call('GMWM_scales_cpp', PACKAGE = 'GMWM', x$nlevels)
  out = structure(list(variance = out[,1],
                         ci_low = out[,2], 
                        ci_high = out[,3], 
                         robust = robust, 
                            eff = eff,
                              p = p,
                         scales = scales), class = "wvar")
  invisible(out)
}

#' @title Print Wavelet Variances
#' @description Displays the summary table of wavelet variance.
#' @method print wvar
#' @param x A \code{wvar} object.
#' @param ... further arguments passed to or from other methods.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(modwt(x))
#' print( out )
print.wvar = function(x, ...){
  mat = matrix(unlist(x[1:3]),ncol=3,byrow=F)
  colnames(mat) = c("Variance", "Low CI", "High CI")
  rownames(mat) = x$scales
  print(mat)
}

#' @title Summary of Wavelet Variances
#' @description Displays the summary table of wavelet variance in addition to CI values and supplied efficiency.
#' @method summary wvar
#' @param object A \code{wvar} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(modwt(x))
#' summary( out )
summary.wvar = function(object, ...){
  name = if(object$robust){
     "robust" 
  }else{
    "classical"
  }
  cat("Results of the wavelet variance calculation using the ",name, " method.\n",sep="")
  if(object$robust){
    cat("Robust was created using efficiency=",object$eff,"\n",sep="")
  }
  
  cat("The confidence interval was generated using (1-",object$p*2,")*100 \n",sep="")
  
  print(object)
}


#' @title Wrapper to ggplot Wavelet Variances Graph
#' @description Creates the wavelet variance graph
#' @method plot wvar
#' @param x A \code{wvar} object.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 graph containing the wavelet variances.
#' @author JJB
#' @seealso \code{\link{autoplot.wvar}}
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(modwt(x))
#' plot( out )
plot.wvar = function(x, ...){
  autoplot.wvar(x)
}

#' @title Graph Wavelet Variances
#' @description Creates the wavelet variance graph
#' @method autoplot wvar
#' @param object A \code{wvar} object.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 graph containing the wavelet variances.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(modwt(x))
#' autoplot( out )
autoplot.wvar = function(object, ...){
  .x=low=high=trans_breaks=trans_format=math_format=NULL
  name = if(object$robust){ "Robust"} else{ "Classic" }
  WV = data.frame(var = object$variance, low = object$ci_low, high = object$ci_high, scale = object$scales)
  CI = ggplot(WV, aes( x = scale, y = low), colour = "#003C7D") + geom_line(linetype = "dotted", colour = "#003C7D") + 
    geom_line(aes(y = high),linetype = "dotted", colour = "#003C7D") +
    geom_line(aes(y = var)) + geom_point(aes(y = var), size = 3) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), alpha = 0.1, colour = "#003C7D") +
    ggtitle(paste0("Haar Wavelet Variance Representation for ", name, " Calculation"))
  
  CI
}




#' @title Compare Wavelet Variances Together
#' @description Creates the wavelet variance graphs with classic and robust together
#' @method autoplot wvarcomp
#' @param object A \code{data.frame} containing both sets of variances.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 graph containing both the wavelet variances.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' Classic = wvar(modwt(x))
#' Robust = wvar(modwt(x), robust=TRUE)
#' compare.WV(Classic = Classic, Robust = Robust, split = FALSE)
autoplot.wvarcomp = function(object, ...){
  
  scales=low1=high1=WV1=low2=high2=WV2=emp=theo=trans_breaks=trans_format=math_format=.x=NULL
  
  WV = data.frame(WV1 = object$WV1, low1 = object$low1, high1 = object$high1, 
                  WV2 = object$WV2, low2 = object$low2, high2 = object$high2, scales = object$scales)
  
  cols = c("LINE1"="#003C7D","LINE2"="#F47F24")
  cols2 = c("LINE2"="#003C7D","LINE1"="#F47F24")
  
  CI = ggplot(WV, aes( x = scales, y = low1)) + 
    geom_line(aes(colour = "LINE1"), linetype = "dotted") +
    geom_line(aes(y = high1, colour = "LINE1"),linetype = "dotted") +
    geom_line(aes(y = WV1, colour = "LINE1")) + geom_point(aes(y = WV1, colour = "LINE1"), size = 3) +
    geom_line(aes(y = low2, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = high2, colour = "LINE2"),linetype = "dotted") +
    geom_line(aes(y = WV2, colour = "LINE2")) + geom_point(aes(y = WV2, colour = "LINE2"), size = 4, shape = 1) +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low1,rev(high1)), x = c(scales,rev(scales))), fill = "#F47F24", alpha = 0.1) +
    geom_polygon(aes(y = c(low2,rev(high2)), x = c(scales,rev(scales))), fill = "#003C7D", alpha = 0.1) +
    ggtitle("Haar Wavelet Variance Representation") +
    theme(legend.key = element_rect(fill=NA), legend.background = element_rect(fill="gray90", size=.5, linetype="dotted"), 
          legend.justification=c(0,0), legend.position=c(0,0.72)) + scale_fill_discrete(name=" ",labels=c(" 95% CI Classical WV  ", " 95% CI Robust WV  ")) +
    scale_colour_manual(name=" ", labels=c("Classical WV Estimator","Robust WV Estimator"), 
                        values = cols2, guide = guide_legend( fill = c("#F47F24","#003C7D") , colour = c("#F47F24","#003C7D"), override.aes = list(shape = c(16,1))))
  
  
  CI
}

#' @title Compare Wavelet Variances on Split
#' @description Creates the wavelet variance graphs with classic and robust 
#' on separate graphs with the same panel
#' @method autoplot wvarcompSplit
#' @param object A \code{data.frame} containing both sets of variances.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing two graphs of the wavelet variance.
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' Classic = wvar(modwt(x))
#' Robust = wvar(modwt(x), robust=TRUE)
#' compare.WV(Classic = Classic, Robust = Robust, split = TRUE)
autoplot.wvarcompSplit = function(object, ...){
  low=high=trans_breaks=trans_format=math_format=.x=NULL
  
  minval = min(c(object$low1,object$low2))
  maxval = max(c(object$high1,object$high2))
  
  WV = data.frame(var = object$WV1, low = object$low1, high = object$high1, scale = object$scales)
  CI1 = ggplot(WV, aes( x = scale, y = low)) + geom_line(linetype = "dotted", colour = "#F47F24") + 
    geom_line(aes(y = high),linetype = "dotted", colour = "#F47F24") +
    geom_line(aes(y = var), colour = "#F47F24") + geom_point(aes(y = var), size = 3, colour = "#F47F24") +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits = c(minval,maxval)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), fill = "#F47F24", alpha = 0.1) +
    ggtitle("Classical WV") 
  
  WV = data.frame(var = object$WV2, low = object$low2, high = object$high2, scale = object$scales)
  CI2 = ggplot(WV, aes( x = scale, y = low), colour = "#003C7D") + geom_line(linetype = "dotted", colour = "#003C7D") + 
    geom_line(aes(y = high),linetype = "dotted", colour = "#003C7D") +
    geom_line(aes(y = var)) + geom_point(aes(y = var), size = 3, colour = "#003C7D") +
    xlab( expression(paste("Scale ", tau))) + ylab( expression(paste("Wavelet variance ", nu))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)), limits = c(minval,maxval)) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_polygon(aes(y = c(low,rev(high)), x = c(scale,rev(scale))), fill = "#003C7D", alpha = 0.1) +
    ggtitle("Robust WV")
  
  multiplot(CI1, CI2, cols=2)
}

#' @title Compare Wavelet Variances
#' @description Compare the estimates given by the classical and robust methods
#'  of calculating the wavelet variance.
#' @usage compare.wvar(Classic, Robust, split = TRUE)
#' @param Classic A \code{wvar} object that has \code{robust=false}.
#' @param Robust A \code{wvar} object that has \code{robust=true}.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE)
#'  or graphed ontop of each other (FALSE)
#' @author JJB
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' Classic = wvar(modwt(x))
#' Robust = wvar(modwt(x), robust=TRUE)
#' # Graph both the classic and the robust estimator on the same graph
#' compare.wvar(Classic = Classic, Robust = Robust, split = FALSE)
#' # Graph both the classic and the robust estimator on split graphs next to each other
#' compare.wvar(Classic = Classic, Robust = Robust, split = TRUE)
compare.wvar = function(Classic, Robust, split = TRUE){
  if(Classic$robust == TRUE){
    if(Robust$robust == TRUE){
      stop("Double robust estimates detected! You must supply different wavelet variance types to compare. e.g. Classic and Robust")
    }else{
      temp = Classic
      Classic = Robust
      Robust = temp
      warning("Object types supplied were reversed (e.g. Robust instead of Classical). We've fixed that for you.")
    }
  }else{
    if(Robust$robust == FALSE){
      stop("Double classic estimates detected! You must supply different wavelet variance types to compare. e.g. Classic and Robust")
    }
  }
  wv = data.frame(WV1 = Classic$variance, low1 = Classic$ci_low, high1 = Classic$ci_high, 
                  WV2 = Robust$variance, low2 = Robust$ci_low, high2 = Robust$ci_high, scales = Classic$scales)
  if (split == TRUE){
    class(wv) = "wvarcompSplit"
  }else{
    class(wv) = "wvarcomp"
  }
  autoplot(wv) 
}
