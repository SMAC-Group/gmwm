#' @title Create a GMWM TS Object
#' @description Setups a time series oriented object that works well with graphing and summary utilities
#' @export
#' @param data A one column matrix, data.frame, or numeric vector.
#' @param name A string that provides an identifier to the data. If not supplied, the name of the first column will be used. Otherwise, it will be registered as "Dataset X"
#' @param freq A number that provides the rate of samples.
#' @param unit A string that contains the unit expression of the frequency. 
#' @return An S3 object with called gts with the following structure:
#' \itemize{
#'  \item{data}{Data set information}
#'  \item{name}{Name of the Dataset}
#'  \item{freq}{Numeric representation of frequency}
#'  \item{Unit}{String representation of the Unit}
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
#' @param object A \code{gts} object
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param line.type A \code{string} that indicates the type of lines.
#' @param line.color A \code{string} that indicates the color of lines.
#' @param point.size A \code{integer} that indicates the size of points on lines.
#' @param point.shape A \code{integer} that indicates the shape of points on lines.
#' @param graph.title A \code{string} that indicates the title of the graph.
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param unit A \code{string} that indicates the unit on x axis.
#' @param title.x.axis A \code{string} that indicates the label on x axis.
#' @param title.y.axis A \code{string} that indicates the label on y axis.
#' @param ... other arguments passed to specific methods
#' @return A ggplot2 panel containing the graph of time series.
#' @author Wenchao
#' @examples
#' x = gen.ts(WN(sigma2=1), 50)
#' x = gts(x)
#' plot(x, background = 'white')
plot.gts = function(object, background = 'grey',  
                    line.type = 'solid', line.color = '#003C7D',
                    point.size = 0, point.shape = 20,
                    graph.title = NULL, graph.title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, unit = object$unit,
                    title.x.axis = 'x', title.y.axis = 'y', ... ){
  
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'grey'
  }
  
  #prepare data
  df = data.frame(y = object$data, x = 1:nrow(object$data))
  
  p = ggplot(data = df, mapping = aes(x = x,y = y)) + geom_line(linetype = line.type, color = line.color) + 
    geom_point(color = line.color, size = point.size, shape = point.shape)
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  if(!is.null(unit)){
    title.x.axis = paste0(title.x.axis, ' (',unit,')')
  }
  
  p = p +  
    xlab(title.x.axis) + ylab(title.y.axis) + ggtitle(graph.title) +
    theme(
      plot.title = element_text(size=graph.title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size))  
  p
}
