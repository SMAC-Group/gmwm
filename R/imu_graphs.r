
#' @title plot Wavelet Variances in separate and combined way
#' @usage plot_imu (object, type, transparence, point.size, color.line, color.point, color.CI, line.type, ...)
#' @param
#'
plot.imu = function(object, type = 'separate', transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ... ){
  
  if (type == 'separate'){
    autoplot.imu6(object, transparence, point.size, color.line, color.point, color.CI, line.type )
  }
  
  if (type == 'combined'){
    autoplot.imu2(object, transparence, point.size, color.line, color.point, color.CI, line.type )
  }
}





#' @title prepare the data for plotting
#' @usage imu2WV(imu)
#'
#'
imu2WV = function(object,...){
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



#' @title plot in separate type: 6 graphs
#'
#'
#'
autoplot.imu6 = function(object, transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ...){
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
  
  quartz(title = "GMWM plot", width = 10, height = 8)
  multiplot(CI)
  
}



#' @title plot in combined way: 2 graphs
#'
#'
#'
#'
autoplot.imu2 = function(object, transparence = 0.1, point.size = 0, color.line = "black", color.point = "#003C7D", color.CI = "#003C7D", line.type = "solid", ...){
  
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

