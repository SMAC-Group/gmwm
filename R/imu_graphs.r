
#' @title Obtain the Wavelet Variance for imu data
#' @description Generates the Wavelet Variance for each column in the data set. 
#' Data set is molten, \code{scales},\code{Axis} and \code{sensor} are used as id variables (see \code{?melt.data.frame})
#' @param object A \code{data.frame} or \code{matrix} that contains 6 columns. 
#' @param gyroscope A \code{vector} that contains the index of column where gyroscope data (Gyro. X, Gyro. Y and Gyro. Z) is placed
#' @param accelerometer A \code{vector} that contains the index of column where accelerometer data (Accel. X, Accel. Y and Accel. Z) is placed
#' @return An object whose class is \code{imu} that contains the formatted data.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' }
imu2WV = function(object, gyroscope = c(1:3), accelerometer = c(4:6)){
  
  if(length(gyroscope) != 3 || length(accelerometer) != 3){
    stop("Gyroscope or Accelerometer parameters must contain exactly 3 column ids each.")
  }
  
  g = gyroscope
  a = accelerometer
  
  # Assume gyroscope columns
  wv1 = wvar(modwt(object[,g[1]]))
  wv2 = wvar(modwt(object[,g[2]]))
  wv3 = wvar(modwt(object[,g[3]]))
  
  # Assume accelerometer columns
  wv4 = wvar(modwt(object[,a[1]]))
  wv5 = wvar(modwt(object[,a[2]]))
  wv6 = wvar(modwt(object[,a[3]]))
  
  temp = data.frame(WV = c(wv1$variance,wv2$variance,wv3$variance,wv4$variance,wv5$variance,wv6$variance),
                      scales = rep(wv1$scales,6),
                      low = c(wv1$ci_low, wv2$ci_low, wv3$ci_low, wv4$ci_low, wv5$ci_low, wv6$ci_low),
                      high = c(wv1$ci_high, wv2$ci_high, wv3$ci_high, wv4$ci_high, wv5$ci_high, wv6$ci_high),
                      Axis = rep(c("X", "Y", "Z", "X", "Y", "Z"), each = length(wv1$variance)),
                      sensor = rep(c("Gyroscope","Accelerometer"), each = 3*length(wv1$variance)))
  
  #melt the data
  imu.df = melt(temp, id.vars = c('scales','Axis','sensor'))
  
  class(imu.df) = "imu"
  return(imu.df)
}

#' @title Plot the Wavelet Variances in separate and combined way
#' @description Creates a graph of the wavelet variance object. The graphs can be either split (features a CI) or combined (all gyroscopes together and all accelerometer).
#' @method plot imu
#' @param x An \code{imu} object (use \code{imu2WV()} to create it)
#' @param separate A \code{boolean} that indicates whether the graphs should be split or combined. 
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param color.line A \code{string} that indicates the color of the line drawn (only work when graphs are seperate)
#' @param color.CI A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.) (only work when graphs are seperate)
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param legend.title A \code{string} that indicates the title of legend (only work when graphs are combined)
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend (only work when graphs are combined)
#' @param legend.title.size An \code{integer} that indicates the size of title on legend (only work when graphs are combined)
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend (only work when graphs are combined)
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
plot.imu = function(x, separate = TRUE, CI = TRUE, transparence = 0.1, color.line = "black", 
                    color.CI = "#003C7D", line.type = "solid", 
                    graph.title = "Haar Wavelet Variance Representation", graph.title.size= 15, 
                    axis.label.size = 13, axis.tick.size = 11, 
                    title.x.axis = expression(paste("Scale ", tau)),
                    title.y.axis = expression(paste("Wavelet Variance ", nu)),
                    legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, 
                    legend.text.size = 13, ... ){
  
  if (separate){
    class(x) = "imu6"
  }
  else{
    class(x) = "imu2"
  }
  
  autoplot(x, CI = CI, transparence = transparence, color.line = color.line, 
           color.CI = color.CI, line.type = line.type, 
           graph.title = graph.title, graph.title.size= graph.title.size, 
           axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
           title.x.axis = title.x.axis,
           title.y.axis = title.y.axis,
           legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
           legend.text.size = legend.text.size)  
}

#' @title Plot in separate type: 6 graphs
#' @description Plot each WV variance in a separate graph
#' @method autoplot imu6
#' @param object An \code{imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param color.line A \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param color.CI A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param ... Additional options
#' @return A panel containing the split graphs of an IMU sensor.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' plot(df)
#' }
autoplot.imu6 = function(object, CI = TRUE, transparence = 0.1, color.line = "black", 
                         color.CI = "#003C7D", line.type = "solid", graph.title = "Haar Wavelet Variance Representation", graph.title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         title.x.axis = expression(paste("Scale ", tau)),
                         title.y.axis = expression(paste("Wavelet Variance ", nu)), ...){
  #require packages: scales
  WV=scales=.x=low=high=NULL
  
  'obj = data.frame(WV = object$WV,
                   scales = object$scales,
                   low = object$low,
                   high = object$high,
                   Axis = object$Axis,
                   sensor = object$sensor)'
  
  #re-construct the data frame
  obj = data.frame(variable = object$variable,
                   value = object$value,
                   scales = object$scales,
                   #low = object$low,
                   #high = object$high,
                   Axis = object$Axis,
                   sensor = object$sensor)

  'p = ggplot(obj, aes(y = variable, x = scales)) +
    geom_line(colour = color.line) +
    geom_point(colour = color.point, linetype = line.type, size = point.size) +
    geom_polygon(aes(y = c(low,rev(high[1:19]),rev(high[1:19 + 19 ]),
                           rev(high[1:19 + 2*19]),rev(high[1:19 + 3*19]),rev(high[1:19 + 4*19]),
                           rev(high[1:19 + 5*19])), x = c(scales,rev(scales))), alpha = transparence , fill = color.CI)'
  
  p = ggplot() +
    geom_line(data = subset(obj, variable == "WV"), mapping = aes(x = scales, y = value), colour = color.line, linetype = line.type) 
    
  if(CI){
    #construct the data frame to plot CI
    low_data_frame = subset(obj, variable == 'low')
    high_value = subset(obj, variable == 'high')$value
    CI_data_frame = data.frame(low_data_frame, high_value)
    
    p = p + geom_ribbon(data = CI_data_frame, mapping = aes(x = scales, ymin = value, ymax = high_value), fill = alpha(color.CI, transparence))
    #color.CI: a hexadecimal color value
  }
  
  p = p + facet_grid(sensor ~ Axis) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(title.x.axis) + ylab(title.y.axis) + ggtitle(graph.title) +
    theme(
      plot.title = element_text(size=graph.title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size))
  
  #multiplot(CI)
  p
}

#' @title Plot in combined type: 2 graphs
#' @description Plot each WV variance in a separate graph
#' @method autoplot imu2
#' @param object An \code{imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
#' @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param graph.title A \code{string} that indicates the title of the graph
#' @param graph.title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark
#' @param title.x.axis A \code{string} that indicates the label on x axis
#' @param title.y.axis A \code{string} that indicates the label on y axis
#' @param legend.title A \code{string} that indicates the title of legend
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend
#' @param ... Additional options
#' @return A panel containing the split graphs of an IMU sensor.
#' @examples
#' \dontrun{
#' data(imu)
#' df = imu2WV(imu)
#' plot(df, separate=FALSE)
#' }
autoplot.imu2 = function(object, CI = T, transparence = 0.1, line.type = "solid",
                         graph.title = "Haar Wavelet Variance Representation", graph.title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         title.x.axis = expression(paste("Scale ", tau)),
                         title.y.axis = expression(paste("Wavelet Variance ", nu)),
                         legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  
  WV=scales=.x=NULL
  
  #rec-construct data frame
  obj = data.frame(variable = object$variable,
                   value = object$value,
                   scales = object$scales,
                   #low = object$low,
                   #high = object$high,
                   Axis = object$Axis,
                   sensor = object$sensor)

  #construct the data frame to plot the line
  line_data_frame = subset(obj, variable == 'WV')
  
  p = ggplot() + geom_line(data = line_data_frame, mapping = aes(x = scales, y = value, color = Axis), linetype = line.type) 
  
  if(CI){
    #construct the data frame to plot the confidence interval
    low_data_frame = subset(obj, variable == 'low')
    high_value = subset(obj, variable == 'high')$value
    CI_data_frame = data.frame(low_data_frame, high_value)
    
    p = p + geom_ribbon(data = CI_data_frame, mapping = aes(x = scales, ymin = value, ymax = high_value, group = Axis, fill = Axis), alpha = transparence)
  }
    
  p = p + facet_grid(sensor~.) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    xlab(title.x.axis) + ylab(title.y.axis) + ggtitle(graph.title) +
    theme(
      plot.title = element_text(size=graph.title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size)) + 
    scale_colour_hue(name = legend.title)
  
  p
    
}
