#' @title Create a GMWM TS Object
#' @description Setups a time series oriented object that works well with graphing and summary utilities
#' @param data A one column matrix, data.frame, or numeric vector.
#' @param name A string that provides an identifier to the data. If not supplied, the name of the first column will be used. Otherwise, it will be registered as "Dataset X"
#' @param freq A number that provides the rate of samples.
#' @param unit A string that contains the unit expression of the frequency. 
#' @return An S3 object with called gts with the following structure:
#' \itemize{
#'  \item{data} {Data set information}
#'  \item{name} {Name of the Dataset}
#'  \item{freq} {Numeric representation of frequency}
#'  \item{Unit} {String representation of the Unit}
#' }
#' @author JJB
#' @examples
#' m = data.frame(rnorm(5))
#' gts(m)
#' 
#' x = gen.ts(WN(sigma2=1), 50)
#' x = gts(x)
gts = function(data, name = "", freq = 1, unit = "sec"){
  
  # Force data.frame to matrix  
  if (is.data.frame(data)){ 
    data = data.matrix(data)
  }
  
  # Check if the data is in matrix form
  if (is.matrix(data)) {
    ndata = nrow(data)
    colnames(data) = NULL
  } else {
    ndata = length(data)
  }
  if (ndata == 0) {
    stop("Not a valid data object! Please supply a data set with one column that is in either a data.frame, matrix,  or numeric object.")
  }
  
  if(!is(freq,"numeric")){
    stop("freq must be numeric")
  }
  
  if(!is(name,"character") || !is(unit,"character")){
    stop("Name or unit is not a character")
  }
  
  x = structure(list(data = data, name = name, freq = freq, unit = unit), class = "gts")
  invisible(x)
}

#' @title Plot Time Series Data
#' @description This function is implemented with ggplot2.
#' @method plot gts
#' @param x A \code{gts} object
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size A \code{integer} that indicates the size of points on lines.
#' @param point.shape A \code{integer} that indicates the shape of points on lines.
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param to.unit A \code{string} indicating the unit which the data is converted to. It can be sec, min, hour, day, month and year.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.ts(WN(sigma2=1), 50)
#' x = gts(x)
#' plot(x)
plot.gts = function(x, background = 'white',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    title = NULL, title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, to.unit = NULL,
                    axis.x.label = 'Time', axis.y.label = deparse(substitute(x)), ... ){
  x1=y1=NULL
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(class(x)!='gts'){stop('x must be a gts object. Use function gts() to create it.')}
  if(x$freq == 0){stop('freq cannot be zero.')}
  
  #prepare data
  len = length(x$data)
  df = data.frame(y1 = x$data, x1 = 0:(len-1) )
  
  if( is.null(x$unit) && is.null(to.unit)==F ){warning('Unit of x is NULL. Conversion cannot be done.')}
  if( is.null(x$unit) == F ){
    if( is.null(to.unit) == F){
      obj = convert.gts(df$x1, x$unit, to.unit)
      if(obj$converted){
        df$x1 = obj$x
        axis.x.label = paste0(axis.x.label, ' (',to.unit,')')
      }
    }
  }
  #deal with freq
  df$x1 = df$x1/x$freq
  
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


# gts.simplify = function(){
#   
# }

#' @title Convert Unit of Time Series Data
#' @param x A \code{vector}.
#' @param from.unit A \code{string} indicating the unit which the data is converted from.
#' @param to.unit A \code{string} indicating the unit which the data is converted to.
#' @return A list with the following structure:
#' \itemize{
#'  \item{x} {Data.}
#'  \item{converted} {A \code{boolean} indicating whether conversion is made.}
#' }
#' @examples
#' x = seq(60, 3600, 60)
#' convert.gts(x, 'sec', 'min')
#' y = 1:10
#' convert.gts(y, 'hour', 'sec')
convert.gts = function(x, from.unit, to.unit){
  #second, min, hour, day, month, year
  unit = c(se = 1, mi = 2, ho = 3, da = 4, mo = 5, ye = 6)
  
  #assume 1 month = 30 days
  ratio = c(60, 60, 24, 30, 12)
  from.unit_1 = substr(from.unit, 1, 2)
  to.unit_1 = substr(to.unit, 1, 2)
  
  #check unit:
  no_convert = F
  if(from.unit_1 == to.unit_1){no_convert = T}
  if(is.na(unit[from.unit_1]) ) {
    message = paste('No such unit: ', from.unit, '. Available units are sec, min, hour, day, month and year. Conversion is terminated.', sep = '')
    warning(message); no_convert = T}
  if(is.na(unit[to.unit_1]) ) {
    message = paste('No such unit: ', to.unit, '. Available units are sec, min, hour, day, month and year. Conversion is terminated.', sep = '')
    warning(message); no_convert = T}
  
  if(!no_convert){
    if(unit[from.unit_1] > unit[to.unit_1]){
      temp = ratio[unit[to.unit_1]: (unit[from.unit_1]-1)]
      multiplier = prod(temp)
      x = x*multiplier
    }else{
      temp = ratio[unit[from.unit_1]: (unit[to.unit_1]-1) ]
      multiplier = prod(temp)
      x = x/multiplier
    }
  }
  obj = list(x = x, converted = !no_convert)  
  return(obj)
}