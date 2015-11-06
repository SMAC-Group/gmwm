# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Create a GMWM TS Object based on data
#' @description Setups a time series oriented object that works well with graphing and summary utilities
#' @param data A one column \code{matrix}, \code{data.frame}, or a numeric \code{vector}.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @return A \code{gts} object with the structure:
#' \itemize{
#'   \item{x:} {A \code{matirx} that contains x-axis values used to plot}
#'   \item{data:} {A \code{matrix} that contains data for combined processes}
#'   \item{freq:} {Numeric representation of frequency}
#'   \item{unit:} {String representation of the unit}
#'   \item{name:} {Name of the dataset}
#' }
#' @author JJB, Wenchao
#' @examples
#' m = data.frame(rnorm(5))
#' gts(m)
#' 
#' x = gen.gts(WN(sigma2=1), 50)$data
#' x = gts(x)
gts = function(data, freq = 1, unit = NULL, name = ""){
  
  # Force data.frame to matrix  
  if (is.data.frame(data)){ 
    data = data.matrix(data)
  }
  
  # Check if the data is in matrix form
  if (is.matrix(data)) {
    ndata = nrow(data)
    colnames(data) = NULL
    
    #check ncol
    ncolumn = ncol(data)
    if(!ncolumn==1){
      stop('data must have one column.')
    }
  } else {
    ndata = length(data)
  }
  if (ndata == 0) {
    stop("Not a valid data object! Please supply a data set with one column that is in either a data.frame, matrix, or numeric object.")
  }
  
  if(!is(freq,"numeric") || length(freq)!=1){
    stop("freq must be numeric")
  }
  if(freq<=0) {stop('freq must be larger than 0.')}
  
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  x = 0:(ndata-1)
  #x = seq(from = 0, to = (ndata-1), length.out = ndata)
  #x = x/freq ###when generate the object, not deal with freq
  
  out = structure(list(
    x = as.matrix(x),
    data = as.matrix(data),
    name = name, 
    freq = freq, 
    unit = unit), class = 'gts' )
  
  invisible(out)
}


#' @title Create a GMWM TS Object based on model
#' @description Create a \code{gts} object based on a supplied time series model.
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} containing the amount of observations for the time series.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @return A \code{gts} object with the structure:
#' \itemize{
#'   \item{x:} {A \code{matirx} that contains x-axis values to plot}
#'   \item{data:} {A \code{matrix} that contains data for combined processes}
#'   \item{freq:} {Numeric representation of frequency}
#'   \item{unit:} {String representation of the unit}
#'   \item{name:} {Name of the Dataset}
#' }
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' gen.gts(model)
gen.gts = function(model, N = 1000, freq = 1, unit = NULL, name = ""){
  
  # Do we have a valid model?
  if(!(is(model, "ts.model") || is(model, "gmwm"))){
    stop("model must be created from a ts.model or gmwm object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  if(is(model,"gmwm")){
    model = model$model.hat
  }
  
  if(!is(freq,"numeric") || length(freq)!=1){
    stop("freq must be numeric")
  }
  if(freq<=0) {stop('freq must be larger than 0.')}
  
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  #x = seq(from = 0, to = N-1, length.out = N)
  #x = x/freq
  x = 0:(N-1)
  
  # Information Required by GMWM:
  desc = model$desc
  
  obj = model$obj.desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(!model$starting){
    theta = model$theta
    out = .Call('gmwm_gen_model', PACKAGE = 'gmwm', N, theta, desc, obj)
  }else{
    stop("Need to supply initial values within the ts.model object.")
  }
  
  out = structure(list(
    x = as.matrix(x),
    data = as.matrix(out),
    name = name, 
    freq = freq, 
    unit = unit), class = 'gts' )
  
  invisible(out)

}

#' @title Plot Time Series Data
#' @description This function is implemented with ggplot2.
#' @method plot gts
#' @export
#' @param x A \code{gts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.gts(WN(sigma2=1), 50)
#' plot(x)
plot.gts = function(x, to.unit = NULL, background = 'white',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    axis.x.label = NULL, axis.y.label = NULL, ... ){
  
  autoplot.gts(object = x, to.unit = to.unit, background = background,  
               line.type = line.type, line.color = line.color,
               point.size = point.size, point.shape = point.shape,
               title = title, title.size= title.size, 
               axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
               axis.x.label = axis.x.label, axis.y.label = axis.y.label )
  
}

#' @title Plot Time Series Data
#' @description This function is implemented with ggplot2.
#' @method autoplot gts
#' @export
#' @keywords internal
#' @param object A \code{gts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.gts(WN(sigma2=1), 50)
#' autoplot(x)
autoplot.gts = function(object, to.unit = NULL, background = 'white',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    axis.x.label = NULL, axis.y.label = NULL, ... ){
  x1=y1=NULL
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(!is(object,"gts")){stop('object must be a gts object. Use function gts() or gen.gts() to create it.')}
  
  #prepare data
  #len = length(x$data)
  #df = data.frame(y1 = x$data, x1 = 0:(len-1) )
  df = data.frame(y1 = object$data, x1 = object$x )
  
  #deal with freq
  df$x1 = df$x1/object$freq
  
  if( is.null(object$unit) && is.null(to.unit)==F ){warning('Unit of object is NULL. Conversion was not done.')}
  if( is.null(object$unit) == F ){
    if( is.null(to.unit) == F){
      obj = unitConversion(df$x1, object$unit, to.unit)
      if(obj$converted){
        df$x1 = obj$x
        message(paste0('Unit of object is converted from ', object$unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  
  p = ggplot(data = df, mapping = aes(x = x1,y = y1)) + geom_line(linetype = line.type, color = line.color) + 
    geom_point(color = line.color, size = point.size, shape = point.shape)
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  p = p +  
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size))  
  p
}


#' @title Convert Unit of Time Series Data
#' @description Manipulate the units of time to different ones
#' @keywords internal
#' @param x A \code{vector} containing the values on x-axis.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit A \code{string} indicating the unit which the data is converted to.
#' @details
#' The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' Make sure \code{from.unit} and \code{to.unit} are not \code{NULL} before it is passed to this function.
#' @return A \code{list} with the following structure:
#' \itemize{
#'  \item{x:} {Data.}
#'  \item{converted:} {A \code{boolean} indicating whether conversion is made.}
#' }
#' @examples
#' x = seq(60, 3600, 60)
#' unitConversion(x, 'sec', 'min')
#' y = 1:10
#' unitConversion(y, 'hour', 'sec')
unitConversion = function(x, from.unit, to.unit){
  #ns, ms, second, min, hour, day, month, year
  unit = c(ns = 1, ms = 2,se = 3, mi = 4, ho = 5, da = 6, mo = 7, ye = 8)
  
  #assume 1 month = 30 days
  ratio = c(1E6, 1E3, 60, 60, 24, 30, 12)
  from.unit.1 = substr(from.unit, 1, 2)
  to.unit.1 = substr(to.unit, 1, 2)
  
  #check unit:
  no.convert = F
  if(from.unit.1 == to.unit.1){no.convert = T}
  if(is.na(unit[from.unit.1]) ) {
    message = paste('No such unit: ', from.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  if(is.na(unit[to.unit.1]) ) {
    message = paste('No such unit: ', to.unit, '. Supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year". Conversion is terminated.', sep = '')
    warning(message); no.convert = T}
  
  if(!no.convert){
    #print out warning when day is convert to month, or month is converted to day.
    conversionRange = unit[from.unit.1] : unit[to.unit.1]
    if(6 %in% conversionRange && 7 %in% conversionRange){
      warning('Unit conversion might be wrong because this function simply assumes 1 month = 30 days.')
    }
  }
  
  if(!no.convert){
    if(unit[from.unit.1] > unit[to.unit.1]){
      temp = ratio[unit[to.unit.1]: (unit[from.unit.1]-1)]
      multiplier = prod(temp)
      x = x*multiplier
    }else{
      temp = ratio[unit[from.unit.1]: (unit[to.unit.1]-1) ]
      multiplier = prod(temp)
      x = x/multiplier
    }
  }
  obj = list(x = x, converted = !no.convert)  
  return(obj)
}