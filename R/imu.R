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

#' @title Create an IMU Object
#' @description Builds an IMU object that provides the program with gyroscope, accelerometer, and axis information per column in the dataset.
#' @param object A \code{vector} which contains data, or a \code{matrix} or \code{data.frame} which contains the data in each column.
#' @param gyroscope A \code{vector} that contains the index of columns where gyroscope data (such as Gyro. X, Gyro. Y and Gyro. Z) is placed.
#' @param accelerometer A \code{vector} that contains the index of columns where accelerometer data (such as Accel. X, Accel. Y and Accel. Z) is placed.
#' @param axis A \code{vector} that indicates the axises, such as 'X', 'Y', 'Z'.
#' @return An \code{imu} object in the following structure:
#' \itemize{
#'   \item{data:} {A \code{matirx} that contains gyroscope and accelerometer data.}
#'   \item{sensor:} {A \code{vector} that indicates whether data contains gyroscope sensor, accelerometer sensor, or both.}
#'   \item{num.sensor:} {A \code{vector} that indicates how many columns of data are for gyroscope sensor and accelerometer sensor.}
#'   \item{axis:} {axis value such as 'X', 'Y', 'Z'.}
#' }
#' @details 
#' \code{object} can be a numeric vector, matrix or data frame.
#' 
#' \code{gyroscope} and \code{accelerometer} cannot be \code{NULL} at the same time, but it will be fine if one of them is \code{NULL}.
#' Also, in order to plot the graph, the length of \code{gyroscope} and \code{accelerometer} are restricted to be equal.
#' 
#' In \code{axis}, duplicate elements are not alowed. If one of parameters between \code{gyroscope} and \code{accelerometer}
#' is \code{NULL}, specify the axis for each column of data. Check example 1 for help. If both of them are not \code{NULL}, specify the
#' \code{axis} only for one parameter (\code{gyroscope} or \code{accelerometer}). Check example 2 for help.
#' 
#' \code{axis} will be automatically generated if there are less than or equal to 3 axises.
#' 
#' @author JJB, Wenchao
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' 
#' # Example 1
#' test1 = imu(imu6, gyroscope = 1:3, accelerometer = NULL, axis = c('X', 'Y', 'Z'))
#' df1 = wvar.imu(test1)
#' plot(df1)
#' 
#' # Example 2
#' test2 = imu(imu6, gyroscope = 1:2, accelerometer = NULL, axis = c('X', 'Y'))
#' df2 = wvar.imu(test2)
#' plot(df2)
#' 
#' # Example 3
#' test3 = imu(imu6, gyroscope = 1:3, accelerometer = 4:6, axis = c('X', 'Y', 'Z'))
#' df3 = wvar.imu(test3)
#' plot(df3)
#' 
#' # Example 4
#' test4 = imu(imu6, gyroscope = 1:2, accelerometer = 4:5, axis = c('X', 'Y'))
#' df4 = wvar.imu(test4)
#' plot(df4)}
imu = function(object, gyroscope = NULL, accelerometer = NULL, axis = NULL){
  
  #check object
  if(is.null(object) || !(is.numeric(object)||is.data.frame(object)||is.matrix(object)) ) {
    stop('Object must a numeric vector, data frame, or matrix.')
  }
  
  if(is.numeric(object)){
    object = as.matrix(object)
  }
  
  if(is.data.frame(object)){
    object = as.matrix(object)
  }
  colnames(object) = NULL
  
  #check gyro and acce
  gyro = gyroscope;
  acce = accelerometer;
  
  if(is.null(gyro) && is.null(acce)){
    stop('At lease one of parameters (gyroscope or accelerometer) must be not NULL.') 
  }else if(!is.null(gyro) && !is.null(acce)){
    if(length(gyro) != length(acce)){
      stop('Object must have equal number of columns for gyroscope and accelerometer sensor.')
    }
  }
  
  #if( (length(gyro) + length(acce)) != ncol(object)){
  #  stop("The number of columns in object doesn't match your supplied indexes for gyroscope and accelerometer.")
  #}
  
  if(!is.whole( c(gyro, acce)) ){
    stop('Paramater gyroscope and accelerometer must be a vector of integers.')
  }
  
  if(any(gyro> ncol(object)) || any(gyro<1)  ){
    stop('Index for gyroscope is out of bound.')
  }
  if(any(acce> ncol(object)) || any(acce<1)){
    stop('Index for accelerometer is out of bound.')
  }

  ##check axis
  ##if we have gyro and acce
  if(!is.null(gyro) && !is.null(acce)){
    index = gyro
  }else{
    index = c(gyro, acce)
  }
  
  #if axis is supplied, make sure it's good
  if(!is.null(axis)){
    #duplicate elements are not allowed
    if(any( count_models(axis)>1) ){
      stop('axis cannot have duplicated elements.')
    }
    
    if(!is.null(gyro) && !is.null(acce)){  
      if(length(axis) != length(index)){stop('When gyroscope and accelerometer are both not NULL, specify the axis only for one sensor.')}
    }else{
      if(length(axis) != length(index)){stop('If only one parameter between gyroscope and accelerometer are not NULL, specify the axis for each column of data.')}
    }
    
  }else{
    #no axis is supplied. try to generate it automatically.
    if(length(index) == 1) {
      axis = 'X'
    }else if (length(index) == 2) {
      axis = c('X','Y')
    }else if (length(index) == 3) {
      axis = c('X', 'Y', 'Z')
    }else{
      stop('axis cannot be automatically generated. Please supply it by specifying "axis = ...".')
    }
  }
  
  # final data manipulation
  # duplicated codes but should be able to save some time for users
  if(is.null(gyro)){
    #gyro is null, but acce is not null
    out = structure(list(data = object[,acce, drop = F],
                         sensor = "Accelerometer",
                         num.sensor = c(length(gyro), length(acce)),
                         axis = axis), class = "imu")
    
  }else if(is.null(acce)){
    #acce is null, but gyro is not null
    out = structure(list(data = object[,gyro, drop = F],
                         sensor = "Gyroscope",
                         num.sensor = c(length(gyro), length(acce)),
                         axis = axis), class = "imu")
  }else{
    object.gyro = object[,gyro, drop = F]
    object.acce = object[,acce, drop = F]
    
    out = structure(list(data = cbind(object.gyro, object.acce),
                         sensor = c("Gyroscope", "Accelerometer" ),
                         num.sensor = c(length(gyro), length(acce)),
                         axis = axis), class = "imu")
  }

  invisible(out)
}

#' @title Wavelet Variance for IMU Object
#' @description Generates the Wavelet Variance for IMU Object. 
#' @method wvar imu
#' @param x An \code{imu} object. 
#' @param alpha	A \code{double} that indicates the \eqn{\left(1-p\right)*\alpha}{(1-p)*alpha} confidence level.
#' @param robust A \code{boolean} that triggers the use of the robust estimate.
#' @param eff A \code{double} that indicates the efficiency as it relates to an MLE.
#' @return A \code{wvar.imu} object which can be plotted directly.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyroscope = 1:3, accelerometer = 4:6)
#' df = wvar.imu(test)
#' }
wvar.imu = function(x, alpha = 0.05, robust = F, eff = 0.6){
  
  if(!is(x, "imu") ){
    stop('This function can only operate on imu object.')
  }

  ncols = sum(x$num.sensor)

  obj.list = vector("list", ncols)
  for(i in 1:ncols){
    obj.list[[i]] = wvar(x$data[, i], alpha = alpha, robust = robust, eff = eff)
  }
  
  ##begin: generate the data frame
  total.len = 0
  each.len = numeric(ncols)
  for (i in 1:ncols){
    each.len[i] = length(obj.list[[i]]$variance)
    total.len = total.len + each.len[i]
  }
  
  #Initialize empty data frame with right number of rows
  obj = data.frame(WV = numeric(total.len),
                   scales = numeric(total.len),
                   low = numeric(total.len),
                   high = numeric(total.len),
                   axis = 'AXIS',
                   sensor = 'SENSOR', stringsAsFactors=FALSE)
  
  if(x$num.sensor[2] == 0){ ## only "Gyroscope"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      d = each.len[i]
      
      obj[t:(t+d-1),] = data.frame(WV = obj.list[[i]]$variance,
                                   scales = obj.list[[i]]$scales,
                                   low = obj.list[[i]]$ci_low,
                                   high = obj.list[[i]]$ci_high,
                                   axis = x$axis[i], 
                                   sensor = "Gyroscope",
                                   stringsAsFactors=FALSE)
      t = t +d
    }
    
  }else if(x$num.sensor[1] == 0){ #only "Accelerometer"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      d = each.len[i]
      
      obj[t:(t+d-1),] = data.frame(WV = obj.list[[i]]$variance,
                                   scales = obj.list[[i]]$scales,
                                   low = obj.list[[i]]$ci_low,
                                   high = obj.list[[i]]$ci_high,
                                   axis = x$axis[i], 
                                   sensor = "Accelerometer",
                                   stringsAsFactors=FALSE)
      t = t +d
    }
    
  }else{ # both "Gyroscope" and "Accelerometer"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      if(i <= length(x$axis)){
        temp.axis = x$axis[i]
        temp.sensor = "Gyroscope"
      }else{ 
        temp.axis = x$axis[i-length(x$axis)]
        temp.sensor = "Accelerometer"
      }
      
      d = each.len[i]
      obj[t:(t+d-1),] = data.frame(WV = obj.list[[i]]$variance,
                                   scales = obj.list[[i]]$scales,
                                   low = obj.list[[i]]$ci_low,
                                   high = obj.list[[i]]$ci_high,
                                   axis = temp.axis, 
                                   sensor = temp.sensor,
                                   stringsAsFactors=FALSE)
      t = t +d
    }
  }

#   temp = data.frame(WV = c(wv1$variance,wv2$variance,wv3$variance,wv4$variance,wv5$variance,wv6$variance),
#                     scales = rep(wv1$scales,6),
#                     low = c(wv1$ci_low, wv2$ci_low, wv3$ci_low, wv4$ci_low, wv5$ci_low, wv6$ci_low),
#                     high = c(wv1$ci_high, wv2$ci_high, wv3$ci_high, wv4$ci_high, wv5$ci_high, wv6$ci_high),
#                     Axis = rep(c("X", "Y", "Z", "X", "Y", "Z"), each = length(wv1$variance)),
#                     sensor = rep(c("Gyroscope","Accelerometer"), each = 3*length(wv1$variance)))
  
  class(obj) = "wvar.imu"
  return(obj)
}



#' @title Wrapper Function to Plot the Wavelet Variances of IMU Object
#' @description Creates a graph of the wavelet variance for imu object.
#' @method plot wvar.imu
#' @export
#' @param x A \code{wvar.imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the confidence interval.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options.
#' @return A panel containing the graph of an IMU sensor.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyroscope = 1:3, accelerometer = 4:6)
#' df = wvar.imu(test)
#' plot(df)
#' plot(df, CI = F)
#' plot(df, CI = T, line.color = c('black', 'black'), title.size = 18)
#' }
plot.wvar.imu = function(x, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                         line.type = NULL, point.size = NULL, point.shape = NULL,
                         CI.color = "#003C7D", title = "Haar Wavelet Variance Representation", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                         facet.label.size = 13, facet.label.background = "#003C7D33",
                         scales = "free_y",...){
  
  autoplot.wvar.imu(x, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                    line.type = line.type, point.size = point.size, point.shape = point.shape,
                    CI.color = CI.color, title = title, title.size= title.size, 
                    axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                    axis.x.label = axis.x.label,
                    axis.y.label = axis.y.label, 
                    facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                    scales = scales)  
}


#' @title Plot the Wavelet Variances of IMU Object
#' @description Creates a graph of the wavelet variance for imu object.
#' @method autoplot wvar.imu
#' @export
#' @keywords internal
#' @param object A \code{wvar.imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the confidence interval.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options.
#' @return A panel containing the graph of an IMU sensor.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyroscope = 1:3, accelerometer = 4:6)
#' df = wvar.imu(test)
#' autoplot(df)
#' autoplot(df, CI = F)
#' autoplot(df, CI = T, line.color = c('black', 'black'), title.size = 18)
#' }
autoplot.wvar.imu = function(object, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             line.type = NULL, point.size = NULL, point.shape = NULL,
                             CI.color = "#003C7D", title = "Haar Wavelet Variance Representation", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                             facet.label.size = 13, facet.label.background = "#003C7D33",
                             scales = "free_y",...){
  
  if(!is(object, 'wvar.imu')){
    stop("This function can only operate on the wvar.imu object. Please use wvar.imu() to create it.")
  }
  
  class(object) = 'imu6'
  
#   if (split){
#     class(x) = "imu6"
#     if(is.null(facet.label.background)){
#       facet.label.background = alpha("#003C7D", transparence)
#     }
#   }else{
#     class(x) = "imu2"
#     if(is.null(facet.label.background)){
#       facet.label.background = alpha("#003C7D", transparence)
#     }
#   }
  
  #call the graphical function
  autoplot(object, CI = CI, background = background, transparence = transparence, line.color = line.color, 
           line.type = line.type, point.size = point.size, point.shape = point.shape,
           CI.color = CI.color, title = title, title.size= title.size, 
           axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
           axis.x.label = axis.x.label,
           axis.y.label = axis.y.label, 
           facet.label.size = facet.label.size, facet.label.background = facet.label.background,
           scales = scales)  
  
}

#' @title Plot imu object in split type: 6 graphs
#' @description Plot each WV variance in a split graph
#' @method autoplot imu6
#' @export
#' @keywords internal
#' @param object A \code{wvar.imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the confidence interval.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options.
#' @return A panel containing the split graphs of an IMU sensor.
autoplot.imu6 = function(object, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                         line.type = NULL, point.size = NULL, point.shape = NULL,
                         CI.color = "#003C7D", title = "Haar Wavelet Variance Representation", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                         facet.label.size = 13, facet.label.background = "#003C7D33",
                         scales = "free_y",...){
  #require packages: scales
  #WV=scales=.x=low=high=NULL
  value = variable = low = high = .x = NULL
  
  #legend.label = NULL,
  #legend.title = '', legend.key.size = 1.3, legend.title.size = 13, 
  #legend.text.size = 13,
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(CI){
    params = c('line.color', 'line.type', 'point.size', 'point.shape', 'CI.color')
    requireLength = c(2, 2, 2, 2, 1)
    default = list(NULL, NULL,  NULL, NULL, "#003C7D")
    nullIsFine = c(rep(T,5))
  }else{
    params = c('line.color', 'line.type', 'point.size', 'point.shape')
    requireLength = c(1, 1, 1, 1)
    default = list(NULL, NULL,  NULL, NULL)
    nullIsFine = c(rep(T,4))
  }
  
  for (i in 1:length(params)){
    one_param = params[i]
    if( length(get(one_param))!=requireLength[i]){
      isNull = is.null(get(one_param))
      if(isNull && nullIsFine[i]){}else{
        warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                      'is supplied.','Default setting is used.'))
      }
      assign(one_param, default[[i]])
    }
  }
  
  if(CI){
    #default setting
    #first WV, then CI
    if(is.null(line.color)){
      line.color = c("#003C7D", "#003C7D")
    }
    if(is.null(line.type)){
      line.type = c('solid', 'dotted')
    }
    if(is.null(point.size)){
      point.size = c(0, 0)
    }
    if(is.null(point.shape)){
      point.shape = c(20, 20)
    }
    
    ##legend.label
    #if(is.null(legend.label)){
    #  legend.label = c(expression(paste("Empirical WV ", hat(nu))), expression(paste("CI(", hat(nu)," , 0.95)" )) )
    #}
    
    #change the length to meet the requirement of ggplot2
    if(length(line.color) == 2){
      line.color = c(line.color, line.color[2])
    }
    
    if(length(line.type) == 2){
      line.type = c(line.type, line.type[2])
    }
     
    if(length(point.size) == 2){
      point.size = c(point.size, point.size[2])
    }
    
    if(length(point.shape) == 2){
      point.shape = c(point.shape, point.shape[2])
    }
    
    #breaks = c('WV','low')
    #legend.fill = c(NA, alpha(CI.color,transparence) )
    #legend.linetype = c(line.type[1],'blank' )
    #legend.pointshape = c(point.shape[1], NA )
    
  }else{
    
    if(is.null(line.color)){
      line.color = c("#003C7D")
    }
    if(is.null(line.type)){
      line.type = c('solid')
    }
    if(is.null(point.size)){
      point.size = c(0)
    }
    if(is.null(point.shape)){
      point.shape = c(20)
    }
    #if(is.null(legend.label)){
    #  legend.label = parse(text = c(expression(paste("Empirical WV ", hat(nu))) ) )
    #}
    
    #breaks = c('WV')
    
  }
  
  #re-construct the data frame
  if(CI){
    obj = data.frame(WV = object$WV,
                     scales = object$scales,
                     low = object$low,
                     high = object$high,
                     axis = object$axis,
                     sensor = object$sensor, stringsAsFactors = F)
    
  }else{
    obj = data.frame(WV = object$WV,
                     scales = object$scales,
                     axis = object$axis,
                     sensor = object$sensor, stringsAsFactors = F)
  }
  melt.obj = melt(obj, id.vars = c('scales', 'axis', 'sensor'))
  
  p = ggplot() +
    geom_line(data = melt.obj, mapping = aes(x = scales, y = value, linetype = variable, color = variable)) +
    geom_point(data = melt.obj, mapping = aes(x = scales, y = value, size = variable, shape = variable, color = variable)) +
    
    scale_linetype_manual(values = c(line.type)) +
    scale_shape_manual(values = c(point.shape))+
    scale_size_manual(values = c(point.size)) +
    scale_color_manual(values = c(line.color))
    
    #scale_linetype_manual(name = legend.title, values = c(line.type),breaks = breaks, labels = legend.label ) +
    #scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label)+
    #scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    #scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label)
  
  if(CI){
    #construct the data frame to plot CI
    obj.CI = data.frame(scales = object$scales,
                        low = object$low,
                        high = object$high,
                        axis = object$axis,
                        sensor = object$sensor, stringsAsFactors = F)
    
    p = p + 
      geom_ribbon(data = obj.CI, mapping = aes(x = scales, ymin = low, ymax = high), fill = alpha(CI.color, transparence), show_guide = F)
     # guides(colour = guide_legend(override.aes = list(fill = legend.fill, linetype = legend.linetype, shape = legend.pointshape)))
     #CI.color: a hexadecimal color value
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  p = p + theme(legend.position='none')+ theme(strip.background = element_rect(fill= facet.label.background) )
  
  p = p + facet_grid(sensor ~ axis, scales = scales) +
    
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
      #legend.key.size = unit(legend.key.size, "cm"),
      #legend.text = element_text(size = legend.text.size),  
      #legend.title = element_text(size = legend.title.size),
      #legend.background = element_rect(fill="transparent"),
      #legend.text.align = 0, 
      strip.text = element_text(size = facet.label.size) )
  
  p
}

# @title Plot in combined type: 2 graphs
# @description Plot each WV variance in a split graph
# @method autoplot imu2
# @param object An \code{imu} object
# @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
# @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
# @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph
# @param line.type A \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
# @param title A \code{string} that indicates the title of the graph
# @param title.size An \code{integer} that indicates the size of title.
# @param axis.label.size An \code{integer} that indicates the size of label
# @param axis.tick.size An \code{integer} that indicates the size of tick mark
# @param axis.x.label A \code{string} that indicates the label on x axis
# @param axis.y.label A \code{string} that indicates the label on y axis
# @param facet.label.size An \code{integer} that indicates the size of facet label
# @param legend.title A \code{string} that indicates the title of legend
# @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend 
# @param legend.title.size An \code{integer} that indicates the size of title on legend
# @param legend.text.size An \code{integer} that indicates the size of key label on legend
# @param facet.label.background A \code{string} that indicates the background color of the facet label
# @param ... Additional options
# @return A panel containing the split graphs of an IMU sensor.
# @examples
# \dontrun{
# data(imu)
# df = imu2WV(imu)
# plot(df, split=FALSE)
# }
# autoplot.imu2 = function(object, CI = T, background = 'grey', transparence = 0.1, line.type = "solid",
#                          title = "Haar Wavelet Variance Representation", title.size= 15, 
#                          axis.label.size = 13, axis.tick.size = 11, 
#                          axis.x.label = expression(paste("Scale ", tau)),
#                          axis.y.label = expression(paste("Wavelet Variance ", nu)),
#                          facet.label.size = 13,
#                          legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, 
#                          legend.text.size = 13, facet.label.background = alpha("#003C7D", transparence), ...){
#   
#   WV=scales=.x=NULL
#   
#   #rec-construct data frame
#   obj = data.frame(WV = object$WV,
#                    #value = object$value,
#                    scales = object$scales,
#                    low = object$low,
#                    high = object$high,
#                    axis = object$axis,
#                    sensor = object$sensor)
#   
#   #construct the data frame to plot the line
#   line_data_frame = subset(obj, variable == 'WV')
#   
#   p = ggplot() + geom_line(data = line_data_frame, mapping = aes(x = scales, y = value, color = axis), linetype = line.type) 
#   
#   if(CI){
#     #construct the data frame to plot the confidence interval
#     low_data_frame = subset(obj, variable == 'low')
#     high_value = subset(obj, variable == 'high')$value
#     CI_data_frame = data.frame(low_data_frame, high_value)
#     
#     p = p + geom_ribbon(data = CI_data_frame, mapping = aes(x = scales, ymin = value, ymax = high_value, group = Axis, fill = Axis), alpha = transparence)
#   }
#   
#   if( background == 'white'){
#     p = p + theme_bw() 
#     p = p + theme(strip.background = element_rect(fill= facet.label.background) )
#   }
#   
#   p = p + facet_grid(sensor~.) + 
#     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                   labels = trans_format("log10", math_format(10^.x))) +
#     scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                   labels = trans_format("log10", math_format(10^.x))) +
#     xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
#     theme(
#       plot.title = element_text(size=title.size),
#       axis.title.y = element_text(size= axis.label.size),
#       axis.text.y  = element_text(size= axis.tick.size),
#       axis.title.x = element_text(size= axis.label.size),
#       axis.text.x  = element_text(size= axis.tick.size),
#       legend.key.size = unit(legend.key.size, "cm"),
#       legend.text = element_text(size = legend.text.size),  
#       legend.title = element_text(size = legend.title.size),
#       strip.text = element_text(size = facet.label.size) ) + 
#     scale_colour_hue(name = legend.title)
#   
#   p
#   
# }


#' @title Wrapper to Automatic Model Selection Results of IMU Object
#' @description Creates a graph of the automatic model selection result containing the empirical and theoretical wavelet variances. 
#' @method plot auto.imu
#' @export
#' @param x A \code{auto.imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the confidence interval.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options.
#' @return A panel containing the automatic model selection results of an IMU sensor.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyroscope = 1:3, accelerometer = 4:6, axis = c('X', 'Y', 'Z'))
#' df = auto.imu(test)
#' plot(df)
#' plot(df, CI = F)
#' plot(df, CI = T, line.color = c('black', 'black', 'blue'), title.size = 18)
#' }
plot.auto.imu = function(x, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                         line.type = NULL, point.size = NULL, point.shape = NULL,
                         CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                         facet.label.size = 13, facet.label.background = "#003C7D33",
                         scales = "free_y",...){
  
  autoplot.auto.imu(x, CI = CI, background = background, transparence = transparence, line.color = line.color, 
  line.type = line.type, point.size = point.size, point.shape = point.shape,
  CI.color = CI.color, title = title, title.size= title.size, 
  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
  axis.x.label = axis.x.label,
  axis.y.label = axis.y.label, 
  facet.label.size = facet.label.size, facet.label.background = facet.label.background,
  scales = scales)
  
  
}


#' @title Automatic Model Selection Results of IMU Object
#' @description Creates a graph of the automatic model selection result containing the empirical and theoretical wavelet variances. 
#' @method autoplot auto.imu
#' @export
#' @keywords internal
#' @param object A \code{auto.imu} object
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the confidence interval.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.)
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.)
#' @param point.size A \code{vector} of \code{integer} that indicates the size of points on lines. 
#' @param point.shape A \code{vector} of \code{integer} that indicates the shape of points on lines.
#' @param CI.color A \code{string} that indicates the color of the confidence interval (e.g. black, red, #003C7D, etc.)
#' @param title A \code{string} that indicates the title of the graph.
#' @param title.size An \code{integer} that indicates the size of title.
#' @param axis.label.size An \code{integer} that indicates the size of label.
#' @param axis.tick.size An \code{integer} that indicates the size of tick mark.
#' @param axis.x.label A \code{string} that indicates the label on x axis.
#' @param axis.y.label A \code{string} that indicates the label on y axis.
#' @param facet.label.size An \code{integer} that indicates the size of facet label.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options.
#' @return A panel containing the automatic model selection results of an IMU sensor.
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyroscope = 1:3, accelerometer = 4:6, axis = c('X', 'Y', 'Z'))
#' df = auto.imu(test)
#' autoplot(df)
#' autoplot(df, CI = F)
#' autoplot(df, CI = T, line.color = c('black', 'black', 'blue'), title.size = 18)
#' }
autoplot.auto.imu = function(object, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             line.type = NULL, point.size = NULL, point.shape = NULL,
                             CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                             facet.label.size = 13, facet.label.background = "#003C7D33",
                             scales = "free_y",...){
  value = variable = low = high = .x = NULL
  
  ###0. param checking
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(CI){
    params = c('line.color', 'line.type', 'point.size', 'point.shape', 'CI.color')
    requireLength = c(3, 3, 3, 3, 1)
    default = list(NULL, NULL,  NULL, NULL, "#003C7D")
    nullIsFine = c(rep(T,5))
  }else{
    params = c('line.color', 'line.type', 'point.size', 'point.shape')
    requireLength = c(2, 2, 2, 2)
    default = list(NULL, NULL,  NULL, NULL)
    nullIsFine = c(rep(T,4))
  }
  
  for (i in 1:length(params)){
    one_param = params[i]
    if( length(get(one_param))!=requireLength[i]){
      isNull = is.null(get(one_param))
      if(isNull && nullIsFine[i]){}else{
        warning(paste('Parameter', one_param, 'requires', requireLength[i],'elements,','but', length(get(one_param)),
                      'is supplied.','Default setting is used.'))
      }
      assign(one_param, default[[i]])
    }
  }
  
  ###1. pre-process the object: auto.imu
  
  if(!is(object, "auto.imu") ){
    stop('This function can only operate on auto.imu object.')
  }
  
  #what is num.sensor and ncols
  num.sensor = object[[1]][[2]]$num.sensor
  ncols = sum(num.sensor)
  
  #what is axis
  if(num.sensor[1] == 0 || num.sensor[2] == 0){##only "Accelerometer"/only "Gyroscope"
    axis = rep(0, ncols)
    for(i in 1:ncols){
      axis[i] = object[[i]][[2]]$axis
    }
    
  }else{#both 
    axis = rep(0, ncols/2)
    
    for(i in 1:(ncols/2)){
      axis[i] = object[[i]][[2]]$axis
    }
  }
  
  #assume 
  obj.list = vector("list", ncols)
  for(i in 1:ncols){
    obj.list[[i]] = object[[i]][[2]]
    
  ######---------------------#######
  #obj.list[[i]]$scales = obj.list[[i]]$scales/100
  }
  
  ##begin: generate the data frame
  total.len = 0
  each.len = numeric(ncols)
  for (i in 1:ncols){
    each.len[i] = length(obj.list[[i]]$wv.empir)
    total.len = total.len + each.len[i]
  }
  
  #Initialize empty data frame with right number of rows
  obj = data.frame(scales = numeric(total.len),
                   emp = numeric(total.len), 
                   low = numeric(total.len),
                   high = numeric(total.len),
                   theo = numeric(total.len),
                   axis = 'AXIS',
                   sensor = 'SENSOR', stringsAsFactors=FALSE)
  
  if(num.sensor[2] == 0){ ## only "Gyroscope"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      d = each.len[i]
      
      obj[t:(t+d-1),] = data.frame(scales = obj.list[[i]]$scales,
                                   emp = obj.list[[i]]$wv.empir,
                                   low = obj.list[[i]]$ci.low,
                                   high = obj.list[[i]]$ci.high,
                                   theo = obj.list[[i]]$theo,
                                   axis = axis[i], 
                                   sensor = "Gyroscope",
                                   stringsAsFactors=FALSE)
      t = t +d
    }
    
  }else if(num.sensor[1] == 0){ #only "Accelerometer"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      d = each.len[i]
      
      obj[t:(t+d-1),] = data.frame(scales = obj.list[[i]]$scales,
                                   emp = obj.list[[i]]$wv.empir,
                                   low = obj.list[[i]]$ci.low,
                                   high = obj.list[[i]]$ci.high,
                                   theo = obj.list[[i]]$theo,
                                   axis = axis[i], 
                                   sensor = "Accelerometer",
                                   stringsAsFactors=FALSE)
      t = t +d
    }
    
  }else{ # both "Gyroscope" and "Accelerometer"
    #put data into data frame
    t = 1
    for (i in 1:ncols){
      if(i <= length(axis)){
        temp.axis = axis[i]
        temp.sensor = "Gyroscope"
      }else{ 
        temp.axis = axis[i-length(axis)]
        temp.sensor = "Accelerometer"
      }
      
      d = each.len[i]
      obj[t:(t+d-1),] = data.frame(scales = obj.list[[i]]$scales,
                                   emp = obj.list[[i]]$wv.empir,
                                   low = obj.list[[i]]$ci.low,
                                   high = obj.list[[i]]$ci.high,
                                   theo = obj.list[[i]]$theo,
                                   axis = temp.axis, 
                                   sensor = temp.sensor,
                                   stringsAsFactors=FALSE)
      t = t +d
    }
  }
  
  #-----END OF DATA PROCESSING-------------------------
  object = obj
  
  ###2. auto-select param
  if(CI){
    #default setting
    #first WV, then CI
    if(is.null(line.color)){
      line.color = c("#003C7D", "#003C7D", "#F47F24")
    }
    if(is.null(line.type)){
      line.type = c('solid', 'dotted', 'solid')
    }
    if(is.null(point.size)){
      point.size = c(0, 0, 0)
    }
    if(is.null(point.shape)){
      point.shape = c(20, 20, 1)
    }
    
    ##legend.label
    #if(is.null(legend.label)){
    #  legend.label = c(expression(paste("Empirical WV ", hat(nu))), expression(paste("CI(", hat(nu)," , 0.95)" )) )
    #}
    
    #change the length to meet the requirement of ggplot2
    if(length(line.color) == 3){
      line.color = c(line.color[1:2], line.color[2:3])
    }
    
    if(length(line.type) == 3){
      line.type = c(line.type[1:2], line.type[2:3])
    }
    
    if(length(point.size) == 3){
      point.size = c(point.size[1:2], point.size[2:3])
    }
    
    if(length(point.shape) == 3){
      point.shape = c(point.shape[1:2], point.shape[2:3])
    }
    
    #breaks = c('WV','low')
    #legend.fill = c(NA, alpha(CI.color,transparence) )
    #legend.linetype = c(line.type[1],'blank' )
    #legend.pointshape = c(point.shape[1], NA )
    
  }else{
    
    if(is.null(line.color)){
      line.color = c("#003C7D", "#F47F24")
    }
    if(is.null(line.type)){
      line.type = c('solid', 'solid')
    }
    if(is.null(point.size)){
      point.size = c(0,0)
    }
    if(is.null(point.shape)){
      point.shape = c(20, 1)
    }
    #if(is.null(legend.label)){
    #  legend.label = parse(text = c(expression(paste("Empirical WV ", hat(nu))) ) )
    #}
    
    #breaks = c('WV')
    
  }
  
  ### 3.reconstruct data frame
  #re-construct the data frame
  if(CI){
    obj = data.frame(scales = object$scales,
                     emp = object$emp,
                     low = object$low,
                     high = object$high,
                     theo = object$theo,
                     axis = object$axis,
                     sensor = object$sensor, stringsAsFactors = F)
    
  }else{
    obj = data.frame(scales = object$scales,
                     emp = object$emp,
                     theo = object$theo,
                     axis = object$axis,
                     sensor = object$sensor, stringsAsFactors = F)
  }
  melt.obj = melt(obj, id.vars = c('scales', 'axis', 'sensor'))
  
  ### 4. start to plot
  p = ggplot() +
    #geom_line(data = melt.obj, mapping = aes(x = scales, y = value, linetype = variable, color = variable), size = 1) +
    geom_line(data = melt.obj, mapping = aes(x = scales, y = value, linetype = variable, color = variable)) +
    geom_point(data = melt.obj, mapping = aes(x = scales, y = value, size = variable, shape = variable, color = variable)) +
    
    scale_linetype_manual(values = c(line.type)) +
    scale_shape_manual(values = c(point.shape))+
    scale_size_manual(values = c(point.size)) +
    scale_color_manual(values = c(line.color))
  
  if(CI){
    #construct the data frame to plot CI
    obj.CI = data.frame(scales = object$scales,
                        low = object$low,
                        high = object$high,
                        axis = object$axis,
                        sensor = object$sensor, stringsAsFactors = F)
    
    p = p + 
      geom_ribbon(data = obj.CI, mapping = aes(x = scales, ymin = low, ymax = high), fill = alpha(CI.color, transparence), show_guide = F)
    # guides(colour = guide_legend(override.aes = list(fill = legend.fill, linetype = legend.linetype, shape = legend.pointshape)))
    #CI.color: a hexadecimal color value
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  p = p + theme(legend.position='none') + theme(strip.background = element_rect(fill= facet.label.background) )
  
  p = p + facet_grid(sensor ~ axis, scales = scales) +
    
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
      #legend.key.size = unit(legend.key.size, "cm"),
      #legend.text = element_text(size = legend.text.size),  
      #legend.title = element_text(size = legend.title.size),
      #legend.background = element_rect(fill="transparent"),
      #legend.text.align = 0, 
      strip.text = element_text(size = facet.label.size) )
  
  p
}


