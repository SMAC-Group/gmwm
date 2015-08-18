#' @title Generate Latent Time Series based on Model
#' @description Create a time series based on a supplied time series model.
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} indicating the amount of observations for the time series.
#' @param x A numeric \code{vector}, 1-column \code{matirx}, or 1-column \code{data frame} containing x-axis values to plot. It must contain \code{N} observations.
#' @param freq A number that provides the rate of samples. Default value is 1.
#' @param unit A string that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A string that provides an identifier to the data. Default value is an empty string.
#' @return A \code{lts} object with the structure:
#' \itemize{
#'   \item{x} {A \code{matirx} that contains x-axis values to plot}
#'   \item{process} {A \code{vector} that contains model names of decomposed and combined processes}
#'   \item{data} {A \code{matrix} that contains data for decomposed and combined processes}
#'   \item{freq} {Numeric representation of frequency}
#'   \item{unit} {String representation of the Unit}
#'   \item{name} {Name of the Dataset}
#' }
#' @author JJB, Wenchao
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' gen.lts(model)
gen.lts = function(model, N = 1000,  x= NULL, freq = 1, unit = NULL, name = ""){
  
  # Do we have a valid model?
  if(!(is(model, "ts.model") || is(model, "gmwm"))){
    stop("model must be created from a ts.model or gmwm object using a supported component (e.g. AR1(), ARMA(p,q), DR(), RW(), QN(), and WN(). ")
  }
  
  if(is(model,"gmwm")){
    model = model$model.hat
  }
  
  if(is.null(x)){x = matrix( 0:(N-1), ncol = 1)}
  
  #check on x: 1) x must be a numeric vector, 1-col matirx or 1-col data frame; 2) it contains N observations
  if(! (is.numeric(x)|| is.data.frame(x) || is.matrix(x))){
    stop("x must be a numeric vector, or 1-column matrix, or 1-column data frame")
  }
  if( (class(x) == "data.frame" || class(x) == "matrix") && ncol(x) != 1) {
    stop("x must be a numeric vector, or 1-column matrix, or 1-column data frame")
  }
  
  if(is.data.frame(x)){
    if(dim(x)[1] != N) {stop(paste0('N = ', N, ', so x must contain ', N, ' observations.') )}
  }else if(length(x) != N){
    stop(paste0('N = ', N, ', so x must contain ', N, ' observations.') )
  }
  
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
    out = .Call('GMWM_gen_lts', PACKAGE = 'GMWM', N, theta, desc, obj)
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

#' @title Plot the Latent Structure in Time Series
#' @description This function is implemented with ggplot2.
#' @method plot lts
#' @param x A \code{lts} object
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param to_unit A \code{string} indicating the unit which the x-axis values is converted to. It can be sec, min, hour, day, month and year.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ncol An \code{integer} that indicates number of columns.
#' @param nrow An \code{integer} that indicates number of rows.
#' @param ... other arguments passed to specific methods.
#' @return A ggplot2 panel containing the graph of time series.
#' @details The unit conversion feature is still under work.
#' @author Wenchao
#' @examples
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' res = gen.lts(model)
#' plot(res)
plot.lts = function(x, to_unit = NULL, background = 'white',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    axis.x.label = 'Time', axis.y.label = object$name, ncol = 1, nrow = NULL, ... ){
  
  autoplot.lts(object = x, to_unit = to_unit, background = background,  
               line.type = line.type, line.color = line.color,
               point.size = point.size, point.shape = point.shape,
               title = title, title.size= title.size, 
               axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
               axis.x.label = axis.x.label, axis.y.label = axis.y.label, ncol = ncol, nrow = nrow)

}

#' @title Plot the Latent Structure in Time Series
#' @description This function is implemented with ggplot2.
#' @method autoplot lts
#' @param object A \code{lts} object
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size An \code{integer} that indicates the size of points on lines.
#' @param point.shape An \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param to_unit A \code{string} indicating the unit which the x-axis values is converted to. It can be sec, min, hour, day, month and year.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ncol An \code{integer} that indicates number of columns
#' @param nrow An \code{integer} that indicates number of rows
#' @param ... other arguments passed to specific methods.
#' @return A ggplot2 panel containing the graph of time series.
#' @details The unit conversion feature is still under work.
#' @author Wenchao
#' @examples
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' res = gen.lts(model)
#' autoplot(res)
autoplot.lts = function(object, to_unit = NULL, background = 'white',  
                        line.type = 'solid', line.color = '#003C7D',
                        point.size = 0, point.shape = 20,
                        title = NULL, title.size= 15, 
                        axis.label.size = 13, axis.tick.size = 11, 
                        axis.x.label = 'Time', axis.y.label = object$name, ncol = 1, nrow = NULL, ...){
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
  
#   if( is.null(x$unit) && is.null(to_unit)==F ){warning('Unit of x is NULL. Conversion cannot be done.')}
#   if( is.null(x$unit) == F ){
#     if( is.null(to_unit) == F){
#       obj = convert.gts(df$x1, x$unit, to_unit)
#       if(obj$converted){
#         df$x1 = obj$x
#         axis.x.label = paste0(axis.x.label, ' (',to_unit,')')
#       }
#     }
#   }
  
  num.desc = length(object$process)
  #A data frame doesn't allow columns to have the same name, but the decomposed processes might be same
  for(i in 1:num.desc) {
    object$process[i] = paste0(object$process[i], paste0(rep(' ',times = (i-1)), collapse = ''))
  }
  
  # prepare data frame to plot
  df = as.data.frame(object$data)
  colnames(df) = object$process
  df$x = object$x
  # melt the data
  melt.df = melt(df, id.vars = c('x'))
  
  p = ggplot(data = melt.df, mapping = aes(x = x,y = value)) + geom_line(linetype = line.type, color = line.color) + 
    geom_point(color = line.color, size = point.size, shape = point.shape)
  
  p = p + facet_wrap(~variable, ncol = ncol, nrow = nrow)
  
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

#' @title Generate a Demo about the Latent Structure in Time Series
#' @description Create a time series based on the supplied model, then generate a demo about its latent structure
#' @param model A \code{ts.model} or \code{gmwm} object containing one of the allowed models.
#' @param N An \code{interger} indicating the amount of observations for the time series.
#' @param x A numeric \code{vector}, 1-column \code{matirx}, or 1-column \code{data frame} containing x-axis values to plot. It must contain \code{N} observations.
#' @param freq A \code{numeric} that provides the rate of samples. Default value is 1.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is an empty string.
#' @author Wenchao
#' @details
#' This function accepts either a \code{ts.model} object (e.g. AR1(phi = .3, sigma2 =1) + WN(sigma2 = 1)) or a \code{gmwm} object.
#' 
#' The unit conversion feature is still under work.
#' @examples
#' # AR
#' set.seed(1336)
#' model = AR1(phi = .99, sigma = 1) + WN(sigma2=1)
#' demo.lts(model)
demo.lts = function(model, N = 1000,  x= NULL, freq = 1, unit = NULL, name = "", ...){
  
  object = gen.lts(model = model, N = N,  x= x, freq = freq, unit = unit, name = name)
  autoplot.lts(object = object, ...)
  
}