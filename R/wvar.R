# Copyright (C) 2014 - 2017  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Wavelet Variance
#' 
#' Calculates the (MODWT) wavelet variance
#' @param x         A \code{vector} with dimensions N x 1, or a \code{lts} object, or a \code{gts} object, or a \code{imu} object. 
#' @param decomp    A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition.
#' @param filter    A \code{string} that specifies what wavelet filter to use. 
#' @param nlevels   An \code{integer} that indicates the level of decomposition. It must be less than or equal to floor(log2(length(x))).
#' @param robust    A \code{boolean} that triggers the use of the robust estimate.
#' @param eff       A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha     A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param freq      A \code{numeric} that provides the rate of samples.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit   A \code{string} indicating the unit which the data is converted to.
#' @param ... Further arguments passed to or from other methods.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @details 
#' If \code{nlevels} is not specified, it is set to \eqn{\left\lfloor {{{\log }_2}\left( {length\left( x \right)} \right)} \right\rfloor}{floor(log2(length(x)))}
#' @author JJB
#' @rdname wvar
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' # Default
#' wvar(x)
#' # Robust
#' wvar(x, robust = TRUE, eff=0.3)
#' # 90% confidence interval
#' wvar(x, alpha = 0.10)
#' 
#' # IMU Object
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyros = 1:3, accels = 4:6, freq = 100)
#' df = wvar.imu(test)
#' }
#' @export
wvar = function(x, ...) {
  UseMethod("wvar")
}

#' @rdname wvar
#' @export
wvar.lts = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  warning('`lts` object is detected. This function can only operate on the combined process.')
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,ncol(x)]

  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.gts = function(x, decomp="modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  freq = attr(x, 'freq')
  unit = attr(x, 'unit')
  x = x[,1]
  
  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.ts = function(x, decomp="modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, to.unit = NULL, ...){
  freq = attr(x, 'tsp')[3]
  unit = NULL

  wvar.default(x, decomp, filter, nlevels, alpha, robust, eff, freq = freq, from.unit = unit, to.unit = to.unit)
}

#' @rdname wvar
#' @export
wvar.default = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = FALSE, eff = 0.6, freq = 1, from.unit = NULL, to.unit = NULL, ...){
  if(is.null(x)){
    stop("`x` must contain a value")
  }else if((is.data.frame(x) || is.matrix(x))){
    if(ncol(x) > 1) stop("There must be only one column of data supplied.")
  }
  
  if(is.null(nlevels)){
    nlevels = floor(log2(length(x)))
  }

  # check freq
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # check unit
  all.units = c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')
  if( (!is.null(from.unit) && !from.unit %in% all.units) || (!is.null(to.unit) && !to.unit %in% all.units) ){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
  }
  
  if(robust) {
    if(eff > 0.99) {
      stop("The efficiency specified is too close to the classical case. Use `robust = FALSE`")
    }
  }
  
  obj =  .Call('gmwm_modwt_wvar_cpp', PACKAGE = 'gmwm',
               signal=x, nlevels=nlevels, robust=robust, eff=eff, alpha=alpha, 
               ci_type="eta3", strWavelet=filter, decomp = decomp)

  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)/freq
  
  # NO unit conversion
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("'from.unit' is NULL. Unit conversion was not done.")
  }
  
  # unit conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      convert.obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (convert.obj$converted) {
        # YES unit conversion
        scales = convert.obj$x
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  if(!is.null(from.unit) && !is.null(to.unit)){ 
    unit = to.unit
  }else{
    unit = from.unit}
  
  create_wvar(obj, decomp, filter, robust, eff, alpha, scales, unit)
}

#' @rdname wvar
#' @export
wvar.imu = function(x, decomp = "modwt", filter = "haar", nlevels = NULL, alpha = 0.05, robust = F, eff = 0.6, to.unit = NULL, ...){

  if(!is.imu(x)){
    stop("`wvar.imu()` requires an IMU Object")
  }
  
  mlevels = floor(log2(nrow(x)))
  
  if(is.null(nlevels)){
    nlevels = mlevels
  }
  if(nlevels > mlevels){
    stop("`nlevels` must be less than ", mlevels,", which is the max number of levels.")
  }
  
  if(robust) {
    if(eff > 0.99) {
      stop("The efficiency specified is too close to the classical case. Use `robust = FALSE`")
    }
  }
  
  # freq conversion
  x.freq = attr(x, 'freq')
  scales = .Call('gmwm_scales_cpp', PACKAGE = 'gmwm', nlevels)/x.freq
  
  # NO unit conversion
  from.unit = attr(x, 'unit')
  if( is.null(from.unit) && is.null(to.unit)==F ){
    warning("The unit of the object is NULL. Unit conversion was not done.")
  }
  
  # unit conversion
  if (!is.null(from.unit)){
    if (!is.null(to.unit)){
      # start_end = c(start, end)
      obj = unitConversion(scales, from.unit = from.unit, to.unit = to.unit)
      
      if (obj$converted) {
        # YES unit conversion
        scales = obj$x
        
        message(paste0('Unit of object is converted from ', from.unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  obj.mat = .Call('gmwm_batch_modwt_wvar_cpp', PACKAGE = 'gmwm', 
                   x, nlevels, robust, eff, alpha, ci_type="eta3", strWavelet="haar", decomp)
  
  # The correct unit
  if(!is.null(from.unit) && !is.null(to.unit)){ 
    unit = to.unit
  }else{
    unit = from.unit}

  # put data into the expected format
  obj.list = lapply(obj.mat, FUN = create_wvar, 
                    decomp = decomp, filter = filter, 
                    robust = robust, eff = eff, 
                    alpha = alpha, scales = scales, unit = unit)
  
  sensor = attr(x, 'sensor')
  x.axis = attr(x, 'axis')
  stype = attr(x, 'stype')
  
  out = structure(list(dataobj = obj.list, #most info is stored in obj.list
                       axis = x.axis,
                       sensor = sensor,
                       stype = stype,
                       freq = x.freq), class="wvar.imu")
  
  out
}


#' Create a \code{wvar} object
#' 
#' Structures elements into a \code{wvar} object
#' @param obj    A \code{matrix} with dimensions N x 3, that contains the wavelet variance, low ci, hi ci.
#' @param decomp A \code{string} that indicates whether to use the "dwt" or "modwt" decomposition
#' @param filter A \code{string} that specifies the type of wavelet filter used in the decomposition
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff    A \code{double} that indicates the efficiency as it relates to an MLE.
#' @param alpha  A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level 
#' @param scales A \code{vec} that contains the amount of decomposition done at each level.
#' @param unit   A \code{string} that contains the unit expression of the frequency.
#' @return A \code{list} with the structure:
#' \describe{
#'   \item{"variance"}{Wavelet Variance}
#'   \item{"ci_low"}{Lower CI}
#'   \item{"ci_high"}{Upper CI}
#'   \item{"robust"}{Robust active}
#'   \item{"eff"}{Efficiency level for Robust}
#'   \item{"alpha"}{p value used for CI}
#'   \item{"unit"}{String representation of the unit}
#' }
#' @keywords internal
create_wvar = function(obj, decomp, filter, robust, eff, alpha, scales, unit){
  structure(list(variance = obj[,1],
                       ci_low = obj[,2], 
                       ci_high = obj[,3], 
                       robust = robust, 
                       eff = eff,
                       alpha = alpha,
                       scales = scales,
                       decomp = decomp,
                       unit = unit,
                       filter = filter), class = "wvar")
}

#' Print Wavelet Variances
#' 
#' Displays the summary table of wavelet variance.
#' @method print wvar
#' @export
#' @keywords internal
#' @param x A \code{wvar} object.
#' @param ... further arguments passed to or from other methods.
#' @author JJB
#' @return Summary table
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = wvar(x)
#' print( out )
print.wvar = function(x, ...){
  mat = matrix(unlist(x[1:3]),ncol=3,byrow=F)
  colnames(mat) = c("Variance", "Low CI", "High CI")
  rownames(mat) = x$scales
  print(mat)
}


#' Print Wavelet Variances for \code{imu} Object
#' 
#' Displays the summary table of wavelet variance for \code{imu} Object.
#' @method print wvar.imu
#' @export
#' @keywords internal
#' @param x A \code{wvar.imu} object.
#' @param ... further arguments passed to or from other methods.
#' @return Summary table
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#' install_imudata()
#' library("imudata")
#' }
#' 
#' data(imu6)
#' 
#' test = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('x', 'y', 'z'), freq = 100)
#' x = wvar(test)
#' print( x )
#' }
print.wvar.imu = function(x, ...){
  dataobj = x$dataobj
  n.obj = length(dataobj)
  
  cat('Level of Decomposition: ', length(x$dataobj[[1]]$scales),", Number of Signals: ", length(x$axis), '\n', sep = '')
  
  if(!is.null(x$stype)){
    cat("Sensor:", x$stype,"@", x$freq,"Hz\n")
    
  }else{
    cat("Freq:", x$freq,"Hz\n")
  }
  
  sensor = x$sensor
  axis = x$axis
  comb = paste(sensor, axis)
  
  for(i in 1:n.obj){
    cat(comb[i], ":\n", sep = "")
    print.wvar(dataobj[[i]])
    
    if(i!=n.obj){cat('\n')}
  }
}


#' Summary of Wavelet Variances
#' 
#' Displays the summary table of wavelet variance in addition to CI values and
#' supplied efficiency.
#' @method summary wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @author JJB
#' @examples
#' set.seed(999)
#' x = rnorm(100)
#' out = wvar(x)
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
  
  cat("The confidence interval was generated using (1-",object$alpha,")*100 \n",sep="")
  
  print(object)
}

#' Summary of Wavelet Variances for \code{imu} Object
#' 
#' Displays the summary table of wavelet variance in addition to CI values and
#' supplied efficiency for \code{imu} Object.
#' @method summary wvar.imu
#' @export
#' @keywords internal
#' @param object A \code{wvar.imu} object.
#' @return Summary table and other properties of the object.
#' @param ... additional arguments affecting the summary produced.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#' install_imudata()
#' library("imudata")
#' }
#' 
#' data(imu6)
#' 
#' test = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('x', 'y', 'z'), freq = 100)
#' x = wvar(test)
#' summary(x)
#' }
summary.wvar.imu = function(object, ...){
  name = if(object$dataobj[[1]]$robust){
    "robust" 
  }else{
    "classical"
  }
  cat("Results of the wavelet variance calculation using the ",name, " method.\n",sep="")
  if(object$dataobj[[1]]$robust){
    cat("Robust was created using efficiency=",object$dataobj[[1]]$eff,"\n",sep="")
  }
  
  cat("The confidence interval was generated using (1-",object$dataobj[[1]]$alpha,")*100 \n",sep="")
  
  print(object)
}


#' Wrapper to ggplot Wavelet Variances Graph
#' 
#' Creates the wavelet variance graph
#' @method plot wvar
#' @export
#' @param x A \code{wvar} object.
#' @template CommonParams
#' @return A ggplot2 graph containing the wavelet variances.
#' @note Parameter line.type, line.color, point.size, point.shape, legend.label must contain 2 elements if \code{CI = TRUE}.
#' @author JJB, Wenchao
#' @seealso \code{\link{autoplot.wvar}}
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(x)
#' plot( out )
plot.wvar = function(x, CI = T, transparence = 0.1, background = 'white', bw = F, 
                     CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                     point.size = NULL, point.shape = NULL,
                     title = NA, title.size= 15, 
                     axis.label.size = 13, axis.tick.size = 11, 
                     axis.x.label = expression(paste("Scale ", tau)),
                     axis.y.label = expression(paste("Wavelet Variance ", nu)),
                     legend.title = '',  legend.label = NULL,
                     legend.key.size = 1, legend.title.size = 13, 
                     legend.text.size = 13, ...){
  
  autoplot.wvar(x, CI = CI, transparence = transparence, background = background, bw = bw, 
                CI.color = CI.color, line.type = line.type, line.color = line.color,
                point.size = point.size, point.shape = point.shape,
                title = title, title.size= title.size, 
                axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                axis.x.label = axis.x.label,
                axis.y.label = axis.y.label,
                legend.title = legend.title,  legend.label = legend.label,
                legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                legend.text.size = legend.text.size )
}

#' Graph Wavelet Variances
#' 
#' Creates the wavelet variance graph
#' @method autoplot wvar
#' @export
#' @keywords internal
#' @param object A \code{wvar} object.
#' @template CommonParams
#' @return A ggplot2 graph containing the wavelet variances.
#' @note Parameter line.type, line.color, point.size, point.shape, legend.label must contain 2 elements if \code{CI = TRUE}.
#' @author JJB, Wenchao
#' @examples
#' set.seed(999)
#' x=rnorm(100)
#' out = wvar(x)
#' autoplot( out )
autoplot.wvar = function(object, CI = T, transparence = 0.1, background = 'white', bw = F, 
                         CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                         point.size = NULL, point.shape = NULL,
                         title = NA, title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)),
                         legend.title = '',  legend.label =  NULL,
                         legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  
  #check parameter
  params = 'legend.label'
  if(CI){
    requireLength = 2
    legend.label.default = c(bquote("Empirical WV"~hat(nu)), bquote("CI("*hat(nu)*", "*.(1 - object$alpha)*")" ))
  }else{
    requireLength = 1
    legend.label.default = c(bquote("Empirical WV"~hat(nu)))
  }
  
  default = list(legend.label.default)
  nullIsFine = T
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  
  p = graphingVar(object, CI = CI, transparence = transparence, background = background, bw = bw, 
                CI.color = CI.color, line.type = line.type, line.color = line.color,
                point.size = point.size, point.shape = point.shape,
                title = title, title.size= title.size, 
                axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                axis.x.label = axis.x.label,
                axis.y.label = axis.y.label,
                legend.title = legend.title, legend.label = legend.label,
                legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                legend.text.size = legend.text.size )
  
 
  if (is.na(title)){
    name = if(object$robust){"Robust"} else{ "Classic" }
    p = p +
      ggtitle(paste0("Haar Wavelet Variance Representation for ", name, " Calculation"))
  }
  
  p
  
}

#' Graphical Function for Allan Variance, Wavelet Variance and Hadamard Variance
#' 
#' A generic graphical function for Allan Variance, Wavelet Variance and
#' Hadamard Variance
#' @param object A \code{wvar}, or \code{avar}, or \code{hadam} object.
#' @template CommonParams
#' @keywords internal
#' @return A ggplot2 graphical object.
#' @author Wenchao
graphingVar = function(object, CI = TRUE, transparence = 0.1, background = 'white', bw = FALSE, 
                        CI.color = "#003C7D", line.type = NULL, line.color = NULL,
                        point.size = NULL, point.shape = NULL,
                        title = NA, title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = NULL,
                        axis.y.label = NULL,
                        legend.title = '',  legend.label =  NULL,
                        legend.key.size = 1, legend.title.size = 13, 
                        legend.text.size = 13, ...){
  
  .x=low=high=trans_breaks=trans_format=math_format=value=variable=NULL
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  #check parameter
  params = c('line.type', 'line.color', 'point.size', 'point.shape')
  if(CI){
    requireLength = rep(2, times = 4)
    default = list(c('solid','dotted'), c('#003C7D', '#999999'),  c(5, 0), c(20,46))
  }else{
    requireLength = rep(1, times = 4)
    default = list('solid', '#003C7D',  5, 20)
  }
  nullIsFine = rep(T,4)
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  if(bw){
    if(CI){line.color = c("#000000", "#404040")}else{line.color = c("#000000")}
    CI.color = "grey50"
  }
  
  WV = data.frame(var = object$variance, low = object$ci_low, high = object$ci_high, scale = object$scales)
  
  if(CI){
    #process parameter (insert some values)
    params = params[-5]; from = 2; to = 3; times = 1;
    for(i in 1:length(params)){
      real_param = get(params[i])
      target = real_param[from]
      stuff = rep(target, times)
      one_param = params[i]
      
      assign(one_param, c(real_param, stuff))
    }
    
    #other parameter
    breaks = c('var', 'low')
    legend.color = c(NA, alpha(CI.color, transparence) )
    legend.linetype = c(line.type[1], 'blank')
    legend.pointshape = c(point.shape[1], NA)
    
    # put data in the desired format
    melt.wv = melt(WV, id.vars = 'scale')
    
  }else{
    #other parameter
    breaks = c('var')
    
    # put data in the desired format
    melt.wv = melt(WV, id.vars = 'scale', measure.vars = 'var')
  }
  
  p = ggplot() + geom_line(data = melt.wv, mapping = aes(x = scale, y = value, color = variable, linetype = variable)) +
    geom_point(data = melt.wv, mapping =aes(x = scale, y = value, color = variable, size = variable, shape = variable)) +
    
    scale_linetype_manual(name = legend.title, values = c(line.type), breaks = breaks, labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label)+
    scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  if(CI){
  p = p + geom_ribbon(data = WV, mapping = aes(ymin = low , ymax = high, x = scale, y = NULL), alpha = transparence, fill = CI.color, show.legend = T) +
    guides(colour = guide_legend(override.aes = list(fill = legend.color, linetype = legend.linetype, shape = legend.pointshape)))
  }
  
  if( background == 'white'||bw){
    p = p + theme_bw() 
  }
  
  #decide where to place the legend
  legendPlace = placeLegend(WV$var[1], WV$low[ length(WV$low) ], WV$high[ length(WV$high)])  
  p = p +
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      legend.justification=legendPlace[1:2], legend.position=legendPlace[3:4],
      legend.text.align = 0)  

  p
}

#' Detail Implementation to Compare Wavelet Variances
#' 
#' Compare the estimates given by the classical and robust methods of 
#' calculating the wavelet variance.
#' @export
#' @keywords internal
#' @param object A \code{data frame} that contains data in order to plot
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE)
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines for wavelet variance and the edge of confidence interval, respectively. Length of vector must equal to 2. 
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point
#' @param title A \code{string} that indicates the title of the graph 
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param axis.x.label A \code{string} that indicates the label on x axis
#' @param axis.y.label A \code{string} that indicates the label on y axis
#' @param facet.label.size An \code{integer} that indicates the size of facet label
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend. If not \code{NULL}, length of vector must equal to the number of \code{wvar} objects that are passed in.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param nrow An \code{integer} that indicates number of rows
#' @param ... Additional parameters
#' @author JJB, Wenchao
#' @seealso \code{\link{compare_wvar}}
autoplot.wvarComp = function(object, split = TRUE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             CI.color = NULL, line.type = NULL, point.size = NULL, point.shape = NULL,
                             title = "Haar Wavelet Variance Representation", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)),
                             facet.label.size = 13, facet.label.background = "#003C7D33",
                             legend.label = NULL,
                             legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                             legend.text.size = 13, nrow = 1, ...){
  
  scales=low=high=WV=emp=theo=trans_breaks=trans_format=math_format=.x=dataset=NULL
  
  #p = ggplot(data = obj, mapping = aes(x = scales, y = WV)) + geom_line(mapping = aes(color = dataset), linetype = line.type[1]) + geom_point(mapping = aes(color = dataset), size = point.size, shape = point.shape) +
  p = ggplot(data = object, mapping = aes(x = scales, y = WV)) + geom_line(mapping = aes(color = dataset), linetype = line.type[1]) + geom_point(mapping = aes(color = dataset, size = dataset, shape = dataset) ) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) 
  if (!is.null(line.color)){
    #legend.label should work. Not work here. But it is changed when creating 'obj' (in wrapper function)
    p = p + scale_color_manual(name = legend.title, values = line.color)
  } 
      
  p = p + scale_size_manual(name = legend.title, values = point.size ) +
      scale_shape_manual(name = legend.title, values = point.shape)

  if(CI){
    p = p + 
      geom_line(mapping = aes(y = low, color = dataset), linetype = line.type[2]) + geom_line(mapping = aes(y = high, color = dataset), linetype = line.type[2]) + 
      geom_ribbon(mapping = aes(ymin = low, ymax = high, fill = dataset), alpha = transparence) 
    if(!is.null(CI.color)){
      p = p + scale_fill_manual(name = legend.title, values = alpha(CI.color, transparence))
    }
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  if (split){
    p = p + facet_wrap(~dataset,nrow = nrow) + theme(legend.position="none")
  }else{
    if(is.null(line.color)){
      p = p + scale_colour_hue(name = legend.title)
    }
    if(is.null(CI.color)){
      p = p + scale_fill_discrete(name = legend.title)
    }
  }
  
  p = p +  xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      strip.text = element_text(size = facet.label.size),
      strip.background = element_rect(fill= facet.label.background) )
  
  p
}


#' Compare Wavelet Variances
#' 
#' Compare the estimates given by the classical and robust methods of 
#' calculating the wavelet variance.
#' @param ... Any number of \code{wvar} or \code{wvar.imu} objects.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE)
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param auto.label.wvar A \code{boolean} that indicates whether legend label should indicate the objects are robust or classical.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point.
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point. 
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param axis.x.label A \code{string} that indicates the label on x axis
#' @param axis.y.label A \code{string} that indicates the label on y axis
#' @param units A two-element vector indicating the units of gyroscope and accelerometer sensor. Set it to \code{NULL} if units are not needed. 
#' @param facet.label.size An \code{integer} that indicates the size of facet label
#' @param facet.label.background A \code{string} that indicates the background color of the facet label
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param nrow An \code{integer} that indicates number of rows
#' @author JJB, Wenchao
#' @note
#' If you meet the error "polygon edge not found", RStudio is complaining that you don't have enough 
#' space to plot the graph. You can adjust the graphics window, or open an external window. The 
#' function \code{\link{external_graphs}} can be used and it works for all operating systems.
#' 
#' When \code{wvar.imu} objects are supplied, some parameters, e.g. \code{split} and \code{nrow}, are 
#' invalid, since the graph is plot seperately and put in 2 rows by default.
#' 
#' @examples
#' \dontrun{
#' ## Case1: Supplied objects are \code{wvar}:
#' data1 = gen_gts(1000, AR1(phi = .32, sigma2=.01))
#' data2 = gen_gts(2000, ARMA(ar=c(.8,.1), ma=c(.3), sigma2=1))
#' data3 = gen_gts(4000, AR1(phi = .32, sigma2=1))
#' 
#' wv1 = wvar(data1, robust = TRUE)
#' wv2 = wvar(data2)
#' wv3 = wvar(data3)
#' 
#' compare_wvar(wv1, wv2)
#' compare_wvar(wv1, wv2, CI = FALSE)
#' compare_wvar(wv1, wv2, split = FALSE)
#' compare_wvar(wv1, wv2, wv3, split = FALSE)
#' 
#' # Change default setting
#' color = c('green','red','blue')
#' label = c('1','2','3')
#' compare_wvar(wv1, wv2, wv3, 
#'              line.color = color, CI.color = color,
#'               legend.label = label)
#' compare_wvar(wv1, wv2, wv3, 
#'              line.color = color, CI.color = color, 
#'              legend.label = label, split = FALSE)
#' 
#' ## Case2: Supplied objects are \code{wvar.imu}:
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test1 = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('X', 'Y', 'Z'), freq = 100)
#' wv1 = wvar(test1)
#' 
#' test2 = imu(imu6, gyros = 1, accels = 3:4, axis = c('X','X','Y'), freq = 100)
#' wv2 = wvar(test2, robust = T)
#' 
#' compare_wvar(wv1, wv2)
#' compare_wvar(wv1, wv2, auto.label.wvar = F, legend.label = c('data1', 'data2'))
#' }
compare_wvar = function(..., background = 'white', split = TRUE, CI = TRUE, auto.label.wvar = T, transparence = 0.1, line.color = NULL, 
                        CI.color = NULL, line.type = NULL,  point.size = NULL, point.shape = NULL,
                        title = "Haar Wavelet Variance Representation", title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu)),
                        units = c(bquote(rad^2/s^2), bquote(m^2/s^4)),
                        facet.label.size = 13, facet.label.background = "#003C7D33",
                        legend.label = NULL,
                        legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                        legend.text.size = 13, nrow = 1 ){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Set background to 'white'")
    background = 'white'
  }
  
  obj_list = list(...)
  object.names = as.character(substitute(...()))
  numObj = length(obj_list)
  
  if (numObj == 0){
    stop('At least one wvar object should be given')
    
  }else if (numObj == 1){
    ## just plot
    plot(...)
    
  }else{
    
    # Case1: Supplied objects are 'wvar.imu'
    is.wvar.imu = sapply(obj_list, FUN = is, class2 = 'wvar.imu')
    
    # Case2: Supplied objects are 'wvar'
    is.wvar = sapply(obj_list, FUN = is, class2 = 'wvar')
    
    if(!all(is.wvar == T) && !all(is.wvar.imu == T)  ){
      stop("Supplied objects must be either 'wvar' or 'wvar.imu' objects.")
    }
    
    # check parameter
    if(all(is.wvar.imu == T)){
      # Case1: Supplied objects are 'wvar.imu'
      if(CI){
        params = c('line.color', 'line.type', 'CI.color', 'point.size', 'point.shape', 'legend.label')
        requireLength = c(3*numObj, 3*numObj, numObj, 3*numObj, 3*numObj, numObj)
        default = list(NULL, NULL,  NULL, NULL, NULL, NULL)
        nullIsFine = c(rep(T,6))
        
      }else{
        params = c('line.color', 'line.type', 'point.size', 'point.shape', 'legend.label')
        requireLength = c(numObj, numObj, numObj, numObj, numObj)
        default = list(NULL, NULL, NULL, NULL, NULL)
        nullIsFine = c(rep(T,5))
      }
    }else{
      # Case2: Supplied objects are 'wvar'
      # line.type is set here
      params = c('line.type', 'line.color', 'CI.color', 'point.size', 'point.shape', 'legend.label')
      requireLength = c(2, numObj, numObj, numObj, numObj, numObj)
      default = list(c('solid','dotted'), NULL,  NULL, rep(5, numObj), rep(20, numObj), NULL)
      nullIsFine = c(rep(T,6))
    }
    checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
    
    # legend.label
    if(is.null(legend.label)){
      legend.label = object.names
    }
    
    if(auto.label.wvar){
      
      rob = sapply( obj_list, FUN = function(x){
        if(is(x, 'wvar.imu')){
          control = x$dataobj[[1]]$robust
        }else{control = x$robust }
         
        if(control){'(Robust)'}else{'(Classical)'} }
      )
      
      legend.label = paste(legend.label, rob)
    }
    # make sure legend.label does not have duplicates
    legend.label = addSpaceIfDuplicate(legend.label)
    
    
    if(all(is.wvar.imu == T)){
      
      return( compare_wvar.imu(obj.list = obj_list,
                       background = background, CI = CI, auto.label.wvar = auto.label.wvar, transparence = transparence, 
                       line.color = line.color, CI.color = CI.color, line.type = line.type, point.size = point.size, point.shape = point.shape,
                       title = title, title.size = title.size, 
                       axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                       axis.x.label = axis.x.label,
                       axis.y.label = axis.y.label,
                       units = units,
                       facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                       legend.label = legend.label,
                       legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                       legend.text.size = legend.text.size ) )
      
    }
    
    total.len = 0
    each.len = numeric(numObj)
    for (i in 1:numObj){
      each.len[i] = length(obj_list[[i]]$variance)
      total.len = total.len + each.len[i]
    }
    #Initialize empty data frame with right number of rows
    obj = data.frame(WV = numeric(total.len),
                     scales = numeric(total.len),
                     low = numeric(total.len),
                     high = numeric(total.len),
                     dataset = 'XYZ', stringsAsFactors=FALSE)
    
    #put data into data frame
    t = 1
    for (i in 1:numObj){
      d = each.len[i]
      obj[t:(t+d-1),] = data.frame(WV = obj_list[[i]]$variance,
                                   scales = obj_list[[i]]$scales,
                                   low = obj_list[[i]]$ci_low,
                                   high = obj_list[[i]]$ci_high,
                                   dataset = legend.label[i], stringsAsFactors=FALSE)
      t = t +d
    }
    
    if (numObj == 2 ){
      if(is.null(line.color)){
        line.color = c("#003C7D","#F47F24")
      }
      if(is.null(CI.color)){
        CI.color = c("#003C7D","#F47F24")
      }
    }
    
    autoplot.wvarComp(obj, split = split, CI = CI, background = background, transparence = transparence, line.color =line.color, 
                        CI.color = CI.color, line.type = line.type,  point.size = point.size, point.shape = point.shape,
                        title = title, title.size= title.size, 
                        axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                        axis.x.label = axis.x.label,
                        axis.y.label = axis.y.label,
                        facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                        legend.label = legend.label,
                        legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                        legend.text.size = legend.text.size,
                        nrow = nrow)
  }
  
}


#' Compare Wavelet Variances for \code{imu} Object
#' 
#' Internal Function to compare the wavelet variance for \code{imu} object.
#' @param obj.list A \code{list} of \code{wvar.imu} objects.
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param auto.label.wvar A \code{boolean} that indicates whether legend label should indicate the objects are robust or classical.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param CI.color A \code{vector} of \code{string} that indicates the color of confidence interval. 
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of point.
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of point. 
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param units A two-element vector indicating the units of gyroscope and accelerometer sensor. Set it to \code{NULL} if units are not needed. 
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param legend.title A \code{string} that indicates the title of legend.
#' @param legend.label A \code{vector} of \code{string} that indicates the labels on legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend.
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend.
#' @author Wenchao
#' @keywords internal
#' @seealso \code{\link{compare_wvar}}
compare_wvar.imu = function(obj.list, background = 'white', CI = TRUE, auto.label.wvar = T, transparence = 0.1, line.color = NULL, 
                        CI.color = NULL, line.type = NULL,  point.size = NULL, point.shape = NULL,
                        title = "Haar Wavelet Variance Representation", title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = expression(paste("Scale ", tau)),
                        axis.y.label = expression(paste("Wavelet Variance ", nu)),
                        units = c(bquote(rad^2/s^2), bquote(m^2/s^4)),
                        facet.label.size = 13, facet.label.background = "#003C7D33",
                        legend.label = NULL,
                        legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
                        legend.text.size = 13) {
  scales=comb=low=high=dataset=.x=NULL
  
  n.obj = length(obj.list)
  
  # 1. Check each graphical parameter. Reset it to default setting if user passes wrong values.
  # checking statements are moved to compare_wvar.
  
  # line.color
  if(is.null(line.color)){
    if(n.obj == 2){ line.color = c("#003C7D","#F47F24") }else{ line.color = ggColor(n.obj) }
  }
  # CI.color
  if(CI && is.null(CI.color)){
    CI.color = line.color
  }
  # go back to line.color
  if(CI && length(line.color) == n.obj){
    line.color = rep(line.color, each = 3)
  }
  
  # line.type 
  if(is.null(line.type)){
    if(CI){ 
      line.type = rep( c('solid','dotted', 'dotted'), times = n.obj) 
    }else{
      line.type = rep( 'solid', times = n.obj)
    }
  }
  
  # point.size 
  if(is.null(point.size)){
    if(CI){ 
      point.size = rep( c(0, 0, 0), times = n.obj) 
    }else{
      point.size = rep(0, times = n.obj)}
  }
  
  # point.shape 
  if(is.null(point.shape)){
    if(CI){ 
      point.shape = rep( c(20, 20, 20), times = n.obj) 
    }else{
      point.shape = rep(20, times = n.obj)}
  }
  
  # breaks
  breaks = paste('emp', legend.label)
  
  # 2. Rearrange the data into a data frame which can be passed to next step.
  #how large the data frame should be
  
  total.len = 0
  each.len = vector('list', length = n.obj )
  
  for (i in 1:n.obj ){
    
    each.wvar.len = length(obj.list[[i]]$dataobj[[1]]$variance)
    n.wvar = length(obj.list[[i]]$dataobj)
    each.len[[i]] = rep(each.wvar.len, times = n.wvar) 
    #assume: all 'wvar' in one 'wvar.imu' have same 'scales'
    
  }
  
  total.len = sum( sapply(each.len, FUN = sum) )
  
  #initialize empty data frame with right number of rows
  if (CI) {
    obj = data.frame(scales = numeric(total.len),
                     emp = numeric(total.len), 
                     low = numeric(total.len),
                     high = numeric(total.len),
                     axis = character(total.len),
                     sensor = character(total.len),
                     dataset = character(total.len), stringsAsFactors=FALSE)
  }else{
    obj = data.frame(scales = numeric(total.len),
                     emp = numeric(total.len), 
                     axis = character(total.len),
                     sensor = character(total.len),
                     dataset = character(total.len), stringsAsFactors=FALSE)
  }
  
  t = 1
  for (i in 1:n.obj ){
    
    data.obj = obj.list[[i]]$dataobj
    axis = obj.list[[i]]$axis
    sensor = obj.list[[i]]$sensor
    
    # add units to sensor
    sensor = addUnits(units = units, sensor = sensor)
    
    for( j in 1:length(each.len[[i]])){
      
      d = each.len[[i]][j]
      
      if(CI){
        
        obj[t:(t+d-1),] = data.frame(scales = data.obj[[j]]$scales, 
                                     emp = data.obj[[j]]$variance,
                                     low = data.obj[[j]]$ci_low,
                                     high = data.obj[[j]]$ci_high,
                                     axis = axis[j], 
                                     sensor = sensor[j],
                                     dataset = legend.label[i], stringsAsFactors=FALSE)
      }else{
        
        obj[t:(t+d-1),] = data.frame(scales = data.obj[[j]]$scales, 
                                     emp = data.obj[[j]]$variance,
                                     axis = axis[j], 
                                     sensor = sensor[j],
                                     dataset = legend.label[i], stringsAsFactors=FALSE)
        
      }
      
      t = t +d
    }
  }
  
  
  # 3. Convert data from a wide format to a long format.
  melt.obj = melt(obj, id.vars = c('scales', 'axis', 'sensor', 'dataset'))
  
  # order that line.type/point.size should be supplied
  if(CI){ lines = c('emp', 'low', 'high') } else{ lines = c('emp') }
  combination = expand.grid(lines, legend.label)
  order = paste(combination[,1], combination[,2])
  
  # combine 'dataset' and 'variable' to set line type, point size, etc
  melt.obj = cbind(melt.obj, 
                   comb = 
                     factor( paste(melt.obj$variable, melt.obj$dataset), levels = order) )
  
  # 4. Call graphical functions in ggplot2 to generate the graph.
  p = ggplot() +
    geom_line(data = melt.obj, mapping = aes(x = scales, y = value, linetype = comb, color = comb) ) +
    geom_point(data = melt.obj, mapping = aes(x = scales, y = value, size = comb, shape = comb, color = comb) ) +
    
    scale_linetype_manual(name = legend.title, values = c(line.type), breaks = breaks, labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label) +
    scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label)
  
  if(CI){
    CI.breaks = legend.label
    
    # change the order to plot CI
    obj$dataset = factor(obj$dataset, levels = legend.label)
      
    p = p + 
      geom_ribbon(data = obj, 
                  mapping = aes(x = scales, ymin = low, ymax = high, fill = dataset), alpha = transparence) +
      scale_fill_manual(name = legend.title, values = c(CI.color), breaks = CI.breaks, labels = legend.label) 
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  p = p + facet_grid(sensor ~ axis, scales = 'free_y', labeller = label_parsed) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      legend.background = element_rect(fill="transparent"),
      legend.text.align = 0, 
      strip.text = element_text(size = facet.label.size), 
      strip.background = element_rect(fill= facet.label.background) )
  
  p
  
}

