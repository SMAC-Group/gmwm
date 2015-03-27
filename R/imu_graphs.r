
#' @title Obtain the Wavelet Variance for imu data
#' @description Generates the Wavelet Variance for each column in the data set. 
#' @param object A \code{data.frame} or \code{matrix} that contains 6 columns. 
#' @return A \code{data.frame} that contains the formatted data.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' }
imu2WV = function(object){
  wv1 = wvar(modwt(object[,1]))
  wv2 = wvar(modwt(object[,2]))
  wv3 = wvar(modwt(object[,3]))
  wv4 = wvar(modwt(object[,4]))
  wv5 = wvar(modwt(object[,5]))
  wv6 = wvar(modwt(object[,6]))
  
  imu.df = data.frame(WV = c(wv1$variance,wv2$variance,wv3$variance,wv4$variance,wv5$variance,wv6$variance),
                      scales = rep(wv1$scales,6),
                      low = c(wv1$ci_low, wv2$ci_low, wv3$ci_low, wv4$ci_low, wv5$ci_low, wv6$ci_low),
                      high = c(wv1$ci_high, wv2$ci_high, wv3$ci_high, wv4$ci_high, wv5$ci_high, wv6$ci_high),
                      axis = rep(c("X", "Y", "Z", "X", "Y", "Z"), each = length(wv1$variance)),
                      sensor = rep(c("Gyroscope","Accelerometer"), each = 3*length(wv1$variance)))
  
  class(imu.df) = "imu"
  return(imu.df)
}


#' @title Plot the Wavelet Variances in separate and combined way
#' @description Creates a graph of the wavelet variance object. The graphs can be either split (features a CI) or combined (all gyroscopes together and all accelerometer).
#' @method plot imu
#' @param x An \code{imu} object
#' @param separate A \code{boolean} that indicates whether the graphs should be split or combined. 
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param point.size An \code{integer} that indicates the size of the points.
#' @param color.line A \code{string} that indicates the color of the line drawn.
#' @param color.point A \code{string} that is a hexadecimal color value.
#' @param color.CI A \code{string} that is a hexadecimal color value.
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param ... Additional options
#' @return A panel containing the graph of an IMU sensor.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' plot(df)
#' plot(df, separate=FALSE)
#' }
#'
plot.imu = function(x, separate = TRUE, transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ... ){
  
  if (separate){
    class(x) = "imu6"
  }
  else{
    class(x) = "imu2"
  }
  
  autoplot(x, transparence, point.size, color.line, color.point, color.CI, line.type )  
}

#' @title Plot in separate type: 6 graphs
#' @description Plot each WV variance in a separate graph
#' @method autoplot imu6
#' @param object An \code{imu} object
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param point.size An \code{integer} that indicates the size of the points.
#' @param color.line A \code{string} that indicates the color of the line drawn.
#' @param color.point A \code{string} that is a hexadecimal color value.
#' @param color.CI A \code{string} that is a hexadecimal color value.
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param ... Additional options
#' @return A panel containing the split graphs of an IMU sensor.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' plot(df)
#' }
autoplot.imu6 = function(object, transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ...){
  
  WV=scales=.x=low=high=NULL
  
  obj = data.frame(WV = object$WV,
                   scales = object$scales,
                   low = object$low,
                   high = object$high,
                   axis = object$axis,
                   sensor = object$sensor)
  
  p = ggplot(obj, aes(y = WV, x = scales)) +
    geom_line(colour = color.line) +
    geom_point(colour = color.point, linetype = line.type, size = point.size) +
    geom_polygon(aes(y = c(low,rev(high[1:19]),rev(high[1:19 + 19 ]),
                           rev(high[1:19 + 2*19]),rev(high[1:19 + 3*19]),rev(high[1:19 + 4*19]),
                           rev(high[1:19 + 5*19])), x = c(scales,rev(scales))), alpha = transparence , fill = color.CI)
  
  CI = p + facet_grid(sensor ~ axis) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  multiplot(CI)
  
}

#' @title Plot in combined type: 2 graphs
#' @description Plot each WV variance in a separate graph
#' @method autoplot imu2
#' @param object An \code{imu} object
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param point.size An \code{integer} that indicates the size of the points.
#' @param color.line A \code{string} that indicates the color of the line drawn.
#' @param color.point A \code{string} that is a hexadecimal color value.
#' @param color.CI A \code{string} that is a hexadecimal color value.
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param ... Additional options
#' @return A panel containing the split graphs of an IMU sensor.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' plot(df, separate=FALSE)
#' }
autoplot.imu2 = function(object, transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ...){
  
  WV=scales=.x=NULL
  
  gyro = data.frame(WV = object$WV[1:57],
                    scales = object$scales[1:57],
                    low = object$low[1:57],
                    high = object$high[1:57],
                    axis = object$axis[1:57],
                    sensor = object$sensor[1:57])
  
  accele = data.frame(WV = object$WV[58:114],
                      scales = object$scales[58:114],
                      low = object$low[58:114],
                      high = object$high[58:114],
                      axis = object$axis[58:114],
                      sensor = object$sensor[58:114])
  
  p = ggplot(gyro, aes(y = WV, x = scales))+geom_line( aes(colour = axis))+
    geom_point(colour = color.point, linetype = line.type, size = point.size)
  #geom_polygon(aes(y = c(low[1: length(low)/2],  rev(high[1:19 + 3*19]),rev(high[1:19 + 4*19]),rev(high[1:19 + 5*19])),
  
  CI1 = p + facet_grid(sensor ~ .) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  p = ggplot(accele, aes(y = WV, x = scales))+geom_line( aes(colour = axis))+
    geom_point(colour = color.point, linetype = line.type, size = point.size)
  
  
  CI2 = p + facet_grid(sensor ~ .) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  
  multiplot(CI1, CI2, cols=1)	
  
}

