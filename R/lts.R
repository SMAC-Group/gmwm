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

#' @title Generate Latent Time Series Object Based on Data
#' @description Create a \code{lts} object based on a supplied matrix or data frame.
#' @param data A multiple-column \code{matrix} or \code{data.frame}. It must contain at least 2 columns. The last column must equal to the sum of all previous columns.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @param process A \code{vector} that contains model names of decomposed and combined processes.
#' @return A \code{lts} object with the structure:
#' \itemize{
#'   \item{x:} {A \code{matirx} that contains x-axis values to plot}
#'   \item{process:} {A \code{vector} that contains model names of decomposed and combined processes}
#'   \item{data:} {A \code{matrix} that contains data for decomposed and combined processes}
#'   \item{freq:} {Numeric representation of frequency}
#'   \item{unit:} {String representation of the unit}
#'   \item{name:} {Name of the dataset}
#' }
#' @author Wenchao
#' @examples
#' model1 = AR1(phi = .99, sigma = 1) 
#' model2 = WN(sigma2=1)
#' col1 = gen.gts(model1, N = 1000)$data
#' col2 = gen.gts(model2, N = 1000)$data
#' testMat = cbind(col1, col2, col1+col2)
#' testLts = lts(testMat, unit = 'sec', process = c('AR1', 'WN', 'AR1+WN'))
#' plot(testLts)
lts = function(data, freq = 1, unit = NULL, name = "", process = NULL){
  
  if(!is(data,'matrix') && !is(data,'data.frame')){
    stop('data must be a matrix or data frame.')
  }
  
  # Force data.frame to matrix  
  if (is.data.frame(data)){ 
    data = data.matrix(data)
  }

  #check ncol
  ncolumn = ncol(data)
  if(ncolumn<2){
    stop('data must have at least two columns.')
  }
  
  #check ndata
  ndata = nrow(data)
  if(ndata == 0 ){stop('data contains 0 observation.')}
  
  #check: the last column must equal to the sum of all previous columns
  tolerance = 1E-4
  sumAllPreviousColumns = apply(data[,1:(ncolumn-1),drop = F], MARGIN = 1, sum)
  checkVec = sumAllPreviousColumns - data[,ncolumn]
  if(any(checkVec>tolerance)){
    stop(paste0('The last column of data must equal to the sum of all previous columns. The tolerance is ', tolerance,'.' ))
  }
  
  #check process
  if(!is.null(process)){
    if(length(process) != ncolumn ){
      stop(paste0('data has ', ncolumn, ' processes (including the combined one). You must specify the name of each process in parameter "process".') )
    }
  }else{
    process = c(paste(rep('process', times = ncolumn-1), 1:(ncolumn-1)), 'combined process')
  }
  
  #check freq
  if(!is(freq,"numeric") || length(freq)!=1){
    stop("freq must be numeric")
  }
  if(freq<=0) {stop('freq must be larger than 0.')}
  
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  #x = seq(from = x.lim[1], to = x.lim[2], length.out = ndata)
  x = 0:(ndata-1)
  
  
  out = structure(list(
    x = as.matrix(x),
    process = process,
    data = data,
    name = name, 
    freq = freq, 
    unit = unit), class = 'lts' )
  
  invisible(out)
}


#' @title Generate Latent Time Series Object Based on Model
#' @description Create a \code{lts} object based on a supplied time series model.
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} indicating the amount of observations generated in this function.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @return A \code{lts} object with the structure:
#' \itemize{
#'   \item{x:} {A \code{matirx} that contains x-axis values to plot}
#'   \item{process:} {A \code{vector} that contains model names of decomposed and combined processes}
#'   \item{data:} {A \code{matrix} that contains data for decomposed and combined processes}
#'   \item{freq:} {Numeric representation of frequency}
#'   \item{unit:} {String representation of the Unit}
#'   \item{name:} {Name of the Dataset}
#' }
#' @author JJB, Wenchao
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' gen.lts(model)
gen.lts = function(model, N = 1000, freq = 1, unit = NULL, name = ""){
  
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

  #if(class(x.lim)!="numeric" || length(x.lim)!=2){
  #  stop('x.lim must be a numeric vector with 2 values.')}
  
  #x = seq(from = x.lim[1], to = x.lim[2], length.out = N)
  x = 0:(N-1)
  
  # Information Required by GMWM:
  desc = model$desc
  p = length(desc) # p decomposed processes
  obj = model$obj.desc
  
  # Identifiability issues
  if(any( count_models(desc)[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, the model will have identifiability issues. Please submit a new model.")
  }
  
  if(!model$starting){
    theta = model$theta
    out = .Call('gmwm_gen_lts', PACKAGE = 'gmwm', N, theta, desc, obj)
  }else{
    stop("Need to supply initial values within the ts.model object.")
  }
  
  #name of each process
  comp.desc = c(desc, paste0(desc, collapse = '+'))
  
  out = structure(list(
    x = as.matrix(x),
    process = comp.desc,
    data = out,
    name = name, 
    freq = freq, 
    unit = unit), class = 'lts' )
  
  invisible(out)
}

#' @title Wrapper Function to Plot the Graph of Latent Time Series
#' @description This function is implemented with ggplot2.
#' @method plot lts
#' @export
#' @param x A \code{lts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be "grey" or "white".
#' @param scales Same as \code{scales} in \code{facet_wrap()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free" in this function.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines.
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param ncol An \code{integer} that indicates number of columns.
#' @param nrow An \code{integer} that indicates number of rows.
#' @param ... other arguments passed to specific methods.
#' @return A ggplot2 panel containing the graph of latent time series.
#' @author Wenchao
#' @examples
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' res = gen.lts(model)
#' plot(res)
plot.lts = function(x, to.unit = NULL, background = 'white', scales = 'free',
                    line.type = NULL, line.color = NULL,
                    point.size = NULL, point.shape = NULL,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    axis.x.label = NULL, axis.y.label = NULL, facet.label.size = 13,
                    ncol = 1, nrow = NULL, ... ){
  
  autoplot.lts(object = x, to.unit = to.unit, background = background, scales = scales, 
               line.type = line.type, line.color = line.color,
               point.size = point.size, point.shape = point.shape,
               title = title, title.size= title.size, 
               axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
               axis.x.label = axis.x.label, axis.y.label = axis.y.label, facet.label.size = facet.label.size,
               ncol = ncol, nrow = nrow)

}

#' @title Plot the Latent Time Series Graph
#' @description This function is implemented with ggplot2.
#' @method autoplot lts
#' @export
#' @keywords internal
#' @param object A \code{lts} object
#' @param to.unit A \code{string} indicating the unit which the data is converted to. The supported units are "ns"(nanosecond), "ms"(millisecond), "sec", "min", "hour", "day", "month", and "year".
#' @param background A \code{string} that determines the graph background. It can be "grey" or "white".
#' @param scales Same as \code{scales} in \code{facet_wrap()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free" in this function.
#' @param line.type A \code{vector} of \code{string} that indicates the type of lines.
#' @param line.color A \code{vector} of \code{string} that indicates the color of lines.
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines.
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param ncol An \code{integer} that indicates number of columns.
#' @param nrow An \code{integer} that indicates number of rows.
#' @param ... other arguments passed to specific methods.
#' @return A ggplot2 panel containing the graph of latent time series.
#' @author Wenchao
#' @examples
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' res = gen.lts(model)
#' autoplot(res)
autoplot.lts = function(object, to.unit = NULL, background = 'white', scales = 'free', 
                        line.type = NULL, line.color = NULL,
                        point.size = NULL, point.shape = NULL,
                        title = NULL, title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = NULL, axis.y.label = NULL, facet.label.size = 13,
                        ncol = 1, nrow = NULL, ...){
  x=value=variable=NULL
  #if user wants to specify nrow, then set ncol = NULL
  if( !is.null(nrow) ){
    ncol = NULL
    warning('ncol is set to NULL in case that ncol and nrow have conflict.')}
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(class(object)!='lts'){stop('object must be a lts object. Use function gen.lts() to create it.')}
  
  #check graphical params
  params = c('line.type', 'line.color', 'point.size','point.shape')
  numLabel = length(object$process)
  for(i in 1:length(params)){
    one_param = params[i]
    if( !is.null(get(one_param)) && length(get(one_param))!=numLabel){
      warning(paste('Parameter', one_param, 'requires',numLabel,'elements,','but', length(get(one_param)),
                    'is supplied.','Default setting is used.'))
      assign(one_param,NULL)
    }
  }
  
  if(is.null(line.type)){line.type = rep('solid', numLabel)}
  if(is.null(line.color)){
    #Set1 = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628" ,"#F781BF", "#999999")
    #modulus = numLabel%/% 9
    #remainder = numLabel%% 9
    #line.color = c( rep(Set1, times = modulus), Set1[1:remainder] )
    
    line.color = ggColor(numLabel)
  }
  if(is.null(point.size)){point.size = rep(0, numLabel)}
  if(is.null(point.shape)){point.shape = rep(20, numLabel)}
  
  num.desc = length(object$process)
  #A data frame doesn't allow columns to have the same name, but the decomposed processes might be same
  for(i in 1:num.desc) {
    object$process[i] = paste0(object$process[i], paste0(rep(' ',times = (i-1)), collapse = ''))
  }
  
  # prepare data frame to plot
  df = as.data.frame(object$data)
  colnames(df) = object$process
  #need to deal with freq now
  df$x = object$x/object$freq
  #df$x = object$x
  
  if( is.null(object$unit) && is.null(to.unit)==F ){warning('Unit of object is NULL. Conversion was not done.')}
  if( is.null(object$unit) == F ){
    if( is.null(to.unit) == F){
      obj = unitConversion(df$x, object$unit, to.unit)
      if(obj$converted){
        df$x = obj$x
        message(paste0('Unit of object is converted from ', object$unit, ' to ', to.unit), appendLF = T)
      }
    }
  }
  
  # melt the data
  melt.df = melt(df, id.vars = c('x'))
  
  p = ggplot(data = melt.df, mapping = aes(x = x,y = value)) + geom_line( mapping = aes(linetype = variable, color = variable)) + 
    geom_point(mapping =aes(color = variable, size = variable, shape = variable)) +
    scale_color_manual( values = line.color ) +
    scale_size_manual( values = point.size ) +
    scale_shape_manual( values = point.shape) +
    scale_linetype_manual(values = line.type)
    
  
  p = p + facet_wrap(~variable, ncol = ncol, nrow = nrow, scales = scales)
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  p = p +  
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    theme(
      legend.position="none",
      plot.title = element_text(size= title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      strip.text = element_text(size = facet.label.size))  
  p
  
}

#' @title Generate a Demo about the Latent Time Series
#' @description Create a time series based on the supplied model, then generate a demo about its latent structure
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} indicating the amount of observations for the time series.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @param ... Additional parameters passed to \code{autoplot.lts}
#' @author Wenchao
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' demo.lts(model)
demo.lts = function(model, N = 1000, freq = 1, unit = NULL, name = "", ...){
  
  object = gen.lts(model = model, N = N, freq = freq, unit = unit, name = name)
  autoplot.lts(object = object, ...)
  
}