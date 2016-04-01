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
#' @param data A \code{vector} which contains data, or a \code{matrix} or \code{data.frame} which contains the data in each column.
#' @param gyros A \code{vector} that contains the index of columns where gyroscope data (such as Gyro. X, Gyro. Y and Gyro. Z) is placed.
#' @param accels A \code{vector} that contains the index of columns where accelerometer data (such as Accel. X, Accel. Y and Accel. Z) is placed.
#' @param axis A \code{vector} that indicates the axises, such as 'X', 'Y', 'Z'. Please supply the axises for gyroscope data before that for accelerometer data, if gyroscope data exists.
#' @param freq An \code{integer} that provides the frequency for the data.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is \code{NULL}.
#' @return An \code{imu} object in the following attributes:
#' \describe{
#'   \item{sensor}{A \code{vector} that indicates whether data contains gyroscope sensor, accelerometer sensor, or both.}
#'   \item{num.sensor}{A \code{vector} that indicates how many columns of data are for gyroscope sensor and accelerometer sensor.}
#'   \item{axis}{Axis value such as 'X', 'Y', 'Z'.}
#'   \item{freq}{Observations per second.}
#'   \item{unit}{String representation of the unit.}
#'   \item{name}{Name of the dataset.}
#' }
#' @details 
#' \code{data} can be a numeric vector, matrix or data frame.
#' 
#' \code{gyros} and \code{accels} cannot be \code{NULL} at the same time, but it will be fine if one of them is \code{NULL}.
#' In the new implementation, the length of \code{gyros} and \code{accels} do not need to be equal.
#' 
#' In \code{axis}, duplicate elements are not alowed for each sensor. In the new implementation, please specify the axis for each column of data.
#' \code{axis} will be automatically generated if there are less than or equal to 3 axises for each sensor.
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
#' # Example 1 - Only gyros
#' test1 = imu(imu6, gyros = 1:3, axis = c('X', 'Y', 'Z'), freq = 100)
#' df1 = wvar.imu(test1)
#' plot(df1)
#' 
#' # Example 2 - One gyro and one accelerometer
#' test2 = imu(imu6, gyros = 1, accels = 4, freq = 100)
#' df2 = wvar.imu(test2)
#' plot(df2)
#' 
#' # Example 3 - 3 gyros and 3 accelerometers
#' test3 = imu(imu6, gyros = 1:3, accels = 4:6, axis = 
#'                        c('X', 'Y', 'Z', 'X', 'Y', 'Z'), freq = 100)
#' df3 = wvar.imu(test3)
#' plot(df3)
#' 
#' # Example 4 - Custom axis
#' test4 = imu(imu6, gyros = 1:2, accels = 4:6, axis = 
#'                        c('X', 'Y', 'X', 'Y', 'Z'), freq = 100)
#' df4 = wvar.imu(test4)
#' plot(df4)
#' }
imu = function(data, gyros = NULL, accels = NULL, axis = NULL, freq = NULL, unit = NULL, name = NULL){
  
  # 1. Check object
  if(is.null(data) || !(is.numeric(data)||is.data.frame(data)||is.matrix(data)) ) {
    stop('Data must a numeric vector, data frame, or matrix.')
  }
  
  if(is.numeric(data)){
    data = as.matrix(data)
  }
  
  if(is.data.frame(data)){
    data = as.matrix(data)
  }
  colnames(data) = NULL
  
  # 2. Check gyro and acce
  gyro = gyros
  acce = accels
  
  ngyros = length(gyro)
  nacces = length(acce)
  
  if(is.null(gyro) && is.null(acce)){
    stop("At lease one of parameters ('gyros' or 'accels') must be not NULL.") 
  }
  
  # Merge indices
  index = c(gyro, acce)
  
  if(!is.whole(index)){
    stop("Paramater 'gyros' and 'accels' must be vectors of integers.")
  }
  
  if(any(gyro > ncol(data)) || any(gyro < 1)){
    stop('Index for gyroscope is out of bound.')
  }
  if(any(acce > ncol(data)) || any(acce < 1)){
    stop('Index for accelerometer is out of bound.')
  }
  
  # 3. Check 'axis': if the user supplies the axis, check input to make sure it is 'good'.
  if(!is.null(axis)){
    
    if(length(axis)==((ngyros + nacces)/2) && ngyros!=0 && nacces!=0){
      axis = rep(axis, times = 2)
    }else if (length(axis) != (ngyros + nacces)){
      stop('Please specify the axis for each column of data.')
    }
    
    if (ngyros == 0||nacces == 0){
      if( anyDuplicated(axis) ){
        stop('`axis` cannot have duplicated elements.')
      }
    }else if (anyDuplicated(axis[1:ngyros]) || anyDuplicated(axis[(ngyros+1):length(axis)])){
      stop('For each sensor, `axis` cannot have duplicated elements.')
    }

  }else{
    # if the user doesn't supply the axis, guess number of sensors
    if(ngyros > 0 && nacces > 0){
      naxis = if(ngyros == nacces) ngyros else 0
    }else{
      naxis = if(ngyros != 0) ngyros else nacces
    }
    
    axis = switch(as.character(naxis),
                  '1' = 'X',
                  '2' = c('X','Y'),
                  '3' = c('X','Y','Z'),
                  stop('axis cannot be automatically generated. Please supply it by specifying "axis = ...".')
    )
    
    if(ngyros == nacces){
      axis = rep(axis, times = 2)
    }
    
  }
  
  # 4. Check freq
  if(is.null(freq)){
    freq = 100
    warning("`freq` has not been specified. Setting `imu` data's frequency to 100. \n Please recreate the object if the frequency is incorrect.")
  }
  if(!is(freq,"numeric") || length(freq) != 1){ stop("'freq' must be one numeric number.") }
  if(freq <= 0) { stop("'freq' must be larger than 0.") }
  
  # 5. do not need 'start' and 'end'
  
  # 6. unit = NULL
  if(!is.null(unit)){
    if(!unit %in% c('ns', 'ms', 'sec', 'second', 'min', 'minute', 'hour', 'day', 'mon', 'month', 'year')){
      stop('The supported units are "ns", "ms", "sec", "min", "hour", "day", "month", "year". ')
    }
  }
  
  create_imu(data[,index, drop = F], ngyros, nacces, axis, freq, unit = unit, name = name)
}

#' @title Internal IMU Object Construction
#' @description Internal quick build for imu object.
#' @param data A \code{matrix} with dimensions N x length(index)
#' @param ngyros An \code{integer} containing the number of gyroscopes
#' @param naccess An \code{integer} containing the number of accelerometers
#' @param axis A \code{vector} unique representation of elements e.g. x,y,z or x,y or x.
#' @param freq An \code{integer} that provides the frequency for the data.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is \code{NULL}.
#' @param stype A \code{string} that describes the sensor type. Default value is \code{NULL}.
#' @return An \code{imu} object class.
#' @keywords internal
create_imu = function(data, ngyros, nacces, axis, freq, unit = NULL, name = NULL, stype = NULL){
  
  if(ngyros>0 && nacces>0){
    colnames(data) = paste( c(rep('Gyro.', times = ngyros), rep('Accel.', times = nacces)), axis)
  }else if (ngyros > 0){
    colnames(data) = c(paste(rep('Gyro.', times = ngyros), axis))
  }else{
    colnames(data) = c(paste(rep('Accel.', times = nacces), axis))
  }
   
  out = structure(data, 
                  sensor = c(rep("Gyroscope",ngyros), rep("Accelerometer",nacces)),
                  num.sensor = c(ngyros, nacces),
                  axis = axis,
                  freq = freq,
                  unit = unit,
                  name = name,
                  stype = stype,
                  class = c("imu","matrix"))

}

#' Subset an IMU Object
#' 
#' Enables the IMU object to be subsettable. That is, you can load all the data in and then select certain properties.
#' @export
#' @param x    A \code{imu} object
#' @param i    A \code{integer vector} that specifies the rows to subset. If blank, all rows are selected.
#' @param j    A \code{integer vector} that specifies the columns to subset. Special rules apply see details.
#' @param drop A \code{boolean} indicating whether the structure should be preserved or simplified.
#' @return An \code{imu} object class.
#' @details 
#' When using the subset operator, note that all the Gyroscopes are placed at the front of object 
#' and, then, the Accelerometers are placed.
#' 
#' @examples 
#' \dontrun{
#' if(!require("imudata")){
#' install_imudata()
#' library("imudata")
#' }
#' 
#' data(imu6)
#' 
#' # Create an IMU Object that is full. 
#' ex = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('X', 'Y', 'Z', 'X', 'Y', 'Z'), freq = 100)
#' 
#' # Create an IMU object that has only gyros. 
#' ex.gyro = ex[,1:3]
#' ex.gyro2 = ex[,c("Gyro. X","Gyro. Y","Gyro. Z")]
#' 
#' # Create an IMU object that has only accels. 
#' ex.accel = ex[,4:6]
#' ex.accel2 = ex[,c("Accel. X","Accel. Y","Accel. Z")]
#' 
#' # Create an IMU object with both gyros and accels on axis X and Y
#' ex.b = ex[,c(1,2,4,5)]
#' ex.b2 = ex[,c("Gyro. X","Gyro. Y","Accel. X","Accel. Y")]
#' 
#' }
#' 
'[.imu' = function(x, i, j, drop = FALSE){
  
  axis = attr(x,"axis")
  sensor = attr(x,"sensor")
  num.sensor = attr(x,"num.sensor")
  
  # If j is missing, then it is simply lowering the number of observations!
  if(!missing(j)){
    # Select column names picked by user
    
    if(is(j, "character")){
      nc = j
    }else{
      nc = colnames(x)[j]
    }
    
    # Remove structure to get Gyros/Accels
    g = gsub("\\..*","",nc)
    ng = table(g)
    
    # Remove structure to get at X,Y,Z axis.
    g2 = gsub(".* ","",nc)
    axis = g2
  
    num.sensor = c({if(!is.na(ng["Gyro"])) ng["Gyro"] else 0}, {if(!is.na(ng["Accel"])) ng["Accel"] else 0})
  }
  
  if(drop){
    return(NextMethod("[", drop = TRUE))
  }
  
  create_imu(NextMethod("[", drop = FALSE),
             num.sensor[1], num.sensor[2], axis, attr(x,"freq"), attr(x,"unit"), attr(x,"name"), attr(x,"stype"))
  
}

#' @title Read an IMU Binary File into R
#' 
#' @description 
#' Process binary files within the 
#' 
#' @param file A \code{string} containing file names or paths.
#' @param type A \code{string} that contains a supported IMU type given below.
#' @param unit A \code{string} that contains the unit expression of the frequency. Default value is \code{NULL}.
#' @param name A \code{string} that provides an identifier to the data. Default value is \code{NULL}.
#' @details
#' Currently supports the following IMUs:
#' \itemize{
#' \item IMAR
#' \item LN200
#' \item LN200IG
#' \item IXSEA
#' \item NAVCHIP_INT
#' \item NAVCHIP_FLT
#' }
#' 
#' We hope to soon be able to support delimited files.
#' @return An \code{imu} object that contains 3 gyroscopes and 3 accelerometers in that order.
#' @references
#' Thanks goes to Philipp Clausen of Labo TOPO, EPFL, Switzerland, topo.epfl.ch, Tel:+41(0)21 693 27 55
#' for providing a matlab function that reads in IMUs.
#' This function is a heavily modified port of MATLAB code into Armadillo/C++.
#' @examples
#' \dontrun{
#' # Relative
#' setwd("F:/")
#' 
#' a = read.imu(file = "Documents/James/short_test_data.imu", type = "IXSEA")
#' 
#' # Fixed path
#' b = read.imu(file = "F:/Desktop/short_test_data.imu", type = "IXSEA")
#' }
read.imu = function(file, type, unit = NULL, name = NULL){
  d = .Call('gmwm_read_imu', PACKAGE = 'gmwm', file_path = file, imu_type = type)
  
  create_imu(d[[1]][,-1], 3, 3, c('X','Y','Z','X','Y','Z'), d[[2]][1], unit = unit, name = name, stype = type)
}


#' @title Wrapper Function to Plot the Wavelet Variances of IMU Object
#' @description Creates a graph of the wavelet variance for imu object.
#' @method plot wvar.imu
#' @export
#' @param x A \code{wvar.imu} object.
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE).
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
#' @param legend.title A \code{string} that indicates the title of legend. It is onlly used when \code{split = FALSE}.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. It is onlly used when \code{split = FALSE}.
#' @param legend.title.size An \code{integer} that indicates the size of title on legend. It is onlly used when \code{split = FALSE}.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend. It is onlly used when \code{split = FALSE}.
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
#' test = imu(imu6, gyros = 1:3, accels = 4:6, freq = 100)
#' df = wvar.imu(test)
#' 
#' ## Plot in split way
#' plot(df, split = T)
#' plot(df, split = T, CI = F)
#' plot(df, split = T, CI = T, line.color = c('black', 'black'), title.size = 18)
#' 
#' ## Plot in combined way
#' plot(df, split = F)
#' plot(df, split = F, line.color = c('black', 'green', 'red'), CI.color = c('black', 'green', 'red'))
#' }
plot.wvar.imu = function(x, split = TRUE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                         line.type = NULL, point.size = NULL, point.shape = NULL,
                         CI.color = NULL, title = "Haar Wavelet Variance Representation", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                         facet.label.size = 13, facet.label.background = "#003C7D33",
                         legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, legend.text.size = 13,
                         scales = "free_y",...){
  
  autoplot.wvar.imu(x, split = split, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                    line.type = line.type, point.size = point.size, point.shape = point.shape,
                    CI.color = CI.color, title = title, title.size= title.size, 
                    axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                    axis.x.label = axis.x.label,
                    axis.y.label = axis.y.label, 
                    facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                    legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, legend.text.size = legend.text.size,
                    scales = scales)  
}


#' @title Plot the Wavelet Variances of IMU Object
#' @description Creates a graph of the wavelet variance for imu object.
#' @method autoplot wvar.imu
#' @export
#' @keywords internal
#' @param object A \code{wvar.imu} object
#' @param split A \code{boolean} that indicates whether the graphs should be separate (TRUE) or graphed ontop of each other (FALSE).
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
#' @param legend.title A \code{string} that indicates the title of legend. It is onlly used when \code{split = FALSE}.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. It is onlly used when \code{split = FALSE}.
#' @param legend.title.size An \code{integer} that indicates the size of title on legend. It is onlly used when \code{split = FALSE}.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend. It is onlly used when \code{split = FALSE}.
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
#' test = imu(imu6, gyros = 1:3, accels = 4:6, freq = 100)
#' df = wvar.imu(test)
#' 
#' ## Plot in split way
#' autoplot(df, split = T)
#' autoplot(df, split = T, CI = F)
#' autoplot(df, split = T, CI = T, line.color = c('black', 'black'), title.size = 18)
#' 
#' ## Plot in combined way
#' autoplot(df, split = F)
#' autoplot(df, split = F, line.color = c('black', 'green', 'red'), 
#'          CI.color = c('black', 'green', 'red'))
#' }
autoplot.wvar.imu = function(object, split = TRUE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             line.type = NULL, point.size = NULL, point.shape = NULL,
                             CI.color = NULL, title = "Haar Wavelet Variance Representation", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                             facet.label.size = 13, facet.label.background = "#003C7D33",
                             legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, legend.text.size = 13,
                             scales = "free_y",...){
  
  if(!inherits(object, 'wvar.imu')){
    stop("This function can only operate on the wvar.imu object. Please use wvar.imu() to create it.")
  }
  
  # Handle new wvar.imu data structure
  object=object$plotobj 
  
  if (split){
    if(is.null(CI.color)){
      CI.color = "#003C7D"
    }
    
    #call the graphical function
    autoplot.imu6(object, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                  line.type = line.type, point.size = point.size, point.shape = point.shape,
                  CI.color = CI.color, title = title, title.size= title.size, 
                  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                  axis.x.label = axis.x.label,
                  axis.y.label = axis.y.label, 
                  facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                  scales = scales)  
    
  }else{
    
    #call the graphical function
    autoplot.imu2(object, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                  line.type = line.type, point.size = point.size, point.shape = point.shape,
                  CI.color = CI.color, title = title, title.size= title.size, 
                  axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                  axis.x.label = axis.x.label,
                  axis.y.label = axis.y.label, 
                  facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                  legend.title = legend.title, legend.key.size = legend.key.size, legend.title.size = legend.title.size, legend.text.size = legend.text.size,
                  scales = scales)
    
  }
  
}

#' @title Plot the Wavelet Variances of IMU Object in Split Type
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
  
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
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
    
    # Change the length to meet the requirement of ggplot2
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
  
  if(CI){
    #construct the data frame to plot CI
    obj.CI = data.frame(scales = object$scales,
                        low = object$low,
                        high = object$high,
                        axis = object$axis,
                        sensor = object$sensor, stringsAsFactors = F)
    
    p = p + 
      geom_ribbon(data = obj.CI, mapping = aes(x = scales, ymin = low, ymax = high), fill = alpha(CI.color, transparence), show.legend = F)
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

#' @title Plot the Wavelet Variances of IMU Object in Combined Type
#' @description Plot each WV variance in a combined graph
#' @method autoplot imu2
#' @export
#' @keywords internal
#' @param object An \code{wvar.imu} object.
#' @param CI A \code{boolean} that indicates whether the confidence interval should be plotted.
#' @param background A \code{string} that determines the graph background. It can be \code{'grey'} or \code{'white'}.
#' @param transparence A \code{double} that ranges from 0 to 1 that controls the transparency of the graph.
#' @param line.color A \code{vector} of \code{string} that indicates the color of the line drawn (e.g. black, blue, red, etc.).
#' @param line.type A \code{vector} of \code{string} that indicates the type of line (e.g. solid, dotted, etc.).
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
#' @param legend.title A \code{string} that indicates the title of legend.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend.
#' @param facet.label.background A \code{string} that indicates the background color of the facet label.
#' @param scales Same as \code{scales} in \code{facet_grid()} in \code{ggplot2} package: should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y"). The default is "free_y" in this function.
#' @param ... Additional options
#' @return A panel containing the combined graphs of an IMU sensor.
autoplot.imu2 = function(object, CI = T, background = 'white', transparence = 0.1, 
                         line.color = NULL, line.type = NULL,
                         point.size = NULL, point.shape = NULL, CI.color = NULL,
                         title = "Haar Wavelet Variance Representation", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)),
                         facet.label.size = 13,
                         legend.title = 'Axis', legend.key.size = 1.3, legend.title.size = 13, 
                         legend.text.size = 13, facet.label.background = "#003C7D33", scales = "free_y", ...){
  
  value=low=high=WV=.x=axis=NULL
  
  # S1: Checking statement (Reset it to default setting if user passes wrong values)
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  num.axis = length( unique(object$axis) )
  if(CI){
    params = c('line.color', 'line.type', 'point.size', 'point.shape', 'CI.color')
    requireLength = c(num.axis, num.axis, num.axis, num.axis, num.axis)
    default = list(NULL, NULL,  NULL, NULL, NULL)
    nullIsFine = c(rep(T,5))
  }else{
    params = c('line.color', 'line.type', 'point.size', 'point.shape')
    requireLength = c(num.axis, num.axis, num.axis, num.axis)
    default = list(NULL, NULL,  NULL, NULL)
    nullIsFine = c(rep(T,4))
  }
  
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  # S2: Auto-select parameters, if not provided by users
  if(is.null(line.color)){
    line.color = ggColor(num.axis)
  }
  if(is.null(line.type)){
    line.type = rep('solid', num.axis)
  }
  if(is.null(point.size)){
    point.size = rep(0, num.axis)
  }
  if(is.null(point.shape)){
    point.shape = rep(20, num.axis)
  }
  
  ## Whether give user the right to modify the legend.label: Currently No
  
  if(CI){
    if(is.null(CI.color)){
      CI.color = ggColor(num.axis) # Change line.color will not automatically change CI.color
    }
  }
  
  # S2: Rearrange the data into a data frame which can be passed to next step
  obj = data.frame(WV = object$WV,
                   scales = object$scales,
                   #low = object$low,
                   #high = object$high,
                   axis = object$axis,
                   sensor = object$sensor, stringsAsFactors = F)
  melt.obj = melt(obj, id.vars = c('scales', 'axis', 'sensor'))
  
  # S3: Generate the graph
  p = ggplot() +
    geom_line(data = melt.obj, mapping = aes(x = scales, y = value, linetype = axis, color = axis)) +
    geom_point(data = melt.obj, mapping = aes(x = scales, y = value, size = axis, shape = axis, color = axis)  ) +
    
    scale_linetype_manual(name = legend.title, values = c(line.type)) +
    scale_shape_manual(name = legend.title, values = c(point.shape))+
    scale_size_manual(name = legend.title, values = c(point.size)) +
    scale_color_manual(name = legend.title, values = c(line.color))
  
  if(CI){
    #construct the data frame to plot CI
    obj.CI = data.frame(scales = object$scales,
                        low = object$low,
                        high = object$high,
                        axis = object$axis,
                        sensor = object$sensor, stringsAsFactors = F)
    
    p = p + geom_ribbon(data = obj.CI, mapping = aes(x = scales, ymin = low, ymax = high, 
                                                     group = axis, fill = axis), alpha = transparence) +
      scale_fill_manual(name = legend.title, values = c(CI.color)) 
  }
  
  if( background == 'white'){
    p = p + theme_bw()
  }
  
  p = p + facet_grid(sensor~., scales = scales) + theme(strip.background = element_rect(fill= facet.label.background) ) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    xlab(axis.x.label) + ylab(axis.y.label) + ggtitle(title) +
    
    theme(
      plot.title = element_text(size=title.size),
      axis.title.y = element_text(size= axis.label.size),
      axis.text.y  = element_text(size= axis.tick.size),
      axis.title.x = element_text(size= axis.label.size),
      axis.text.x  = element_text(size= axis.tick.size),
      legend.key.size = unit(legend.key.size, "cm"),
      legend.text = element_text(size = legend.text.size),  
      legend.title = element_text(size = legend.title.size),
      strip.text = element_text(size = facet.label.size) )
  
  p
  
}


#' @title Wrapper to Automatic Model Selection Results of IMU Object
#' @description Creates a graph of the automatic model selection result containing the empirical and theoretical wavelet variances. 
#' @method plot auto.imu
#' @export
#' @param x A \code{auto.imu} object
#' @param process.decomp A \code{boolean} that indicates whether the decomposed processes should be plotted or not.
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
#' @param legend.title A \code{string} that indicates the title of legend. Used only when \code{process.decomp = T}.
#' @param legend.label A \code{character vector} that indicates the labels on legend. Used only when \code{process.decomp = T}.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. Used only when \code{process.decomp = T}.
#' @param legend.title.size An \code{integer} that indicates the size of title on legend. Used only when \code{process.decomp = T}.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend. Used only when \code{process.decomp = T}.
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
#' test = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('X', 'Y', 'Z', 'X', 'Y', 'Z'), freq = 100)
#' df = auto.imu(test)
#' plot(df)
#' plot(df, process.decomp = T)
#' plot(df, CI = F)
#' plot(df, CI = T, line.color = c('black', 'black', 'blue'), title.size = 18)
#' }
plot.auto.imu = function(x, process.decomp = FALSE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                         line.type = NULL, point.size = NULL, point.shape = NULL,
                         CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                         axis.label.size = 13, axis.tick.size = 11, 
                         axis.x.label = expression(paste("Scale ", tau)),
                         axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                         facet.label.size = 13, facet.label.background = "#003C7D33", scales = "free_y",
                         legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                         legend.text.size = 13, ...){
  
  autoplot.auto.imu(x, process.decomp = process.decomp, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                    line.type = line.type, point.size = point.size, point.shape = point.shape,
                    CI.color = CI.color, title = title, title.size= title.size, 
                    axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                    axis.x.label = axis.x.label,
                    axis.y.label = axis.y.label, 
                    facet.label.size = facet.label.size, facet.label.background = facet.label.background, scales = scales,
                    legend.title = legend.title,  legend.label = legend.label, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                    legend.text.size = legend.text.size)
  
  
}


#' @title Automatic Model Selection Results of IMU Object
#' @description Creates a graph of the automatic model selection result containing the empirical and theoretical wavelet variances. 
#' @method autoplot auto.imu
#' @export
#' @keywords internal
#' @param object A \code{auto.imu} object
#' @param process.decomp A \code{boolean} that indicates whether the decomposed processes should be plotted or not.
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
#' @param legend.title A \code{string} that indicates the title of legend. Used only when \code{process.decomp = T}.
#' @param legend.label A \code{character vector} that indicates the labels on legend. Used only when \code{process.decomp = T}.
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. Used only when \code{process.decomp = T}.
#' @param legend.title.size An \code{integer} that indicates the size of title on legend. Used only when \code{process.decomp = T}.
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend. Used only when \code{process.decomp = T}.
#' @param ... Additional options.
#' @return A panel containing the automatic model selection results of an IMU sensor.
#' @author Wenchao
#' @examples
#' \dontrun{
#' if(!require("imudata")){
#'    install_imudata()
#'    library("imudata")
#' }
#' 
#' data(imu6)
#' test = imu(imu6, gyros = 1:3, accels = 4:6, axis = c('X', 'Y', 'Z', 'X', 'Y', 'Z'), freq = 100)
#' df = auto.imu(test)
#' autoplot(df)
#' autoplot(df, process.decomp = T)
#' autoplot(df, CI = F)
#' autoplot(df, CI = T, line.color = c('black', 'black', 'blue'), title.size = 18)
#' }
autoplot.auto.imu = function(object, process.decomp = FALSE, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                              line.type = NULL, point.size = NULL, point.shape = NULL,
                              CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                              axis.label.size = 13, axis.tick.size = 11, 
                              axis.x.label = expression(paste("Scale ", tau)),
                              axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                              facet.label.size = 13, facet.label.background = "#003C7D33",scales = "free_y",
                              legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                              legend.text.size = 13, ...){
  
  ## Common checking shared by two functions
  if( !(background %in% c('grey','gray', 'white')) ){
    warning("Parameter background: No such option. Default setting is used.")
    background = 'white'
  }
  
  if(!is(object, "auto.imu") ){
    stop('This function can only operate on auto.imu object.')
  }
  
  if(!process.decomp){
    
    autoplot.auto.imu1(object, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                      line.type = line.type, point.size = point.size, point.shape = point.shape,
                      CI.color = CI.color, title = title, title.size= title.size, 
                      axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                      axis.x.label = axis.x.label,
                      axis.y.label = axis.y.label, 
                      facet.label.size = facet.label.size, facet.label.background = facet.label.background,
                      scales = scales)
  }else{
    
    #process.decomp == T
    autoplot.auto.imu2(object, CI = CI, background = background, transparence = transparence, line.color = line.color, 
                       line.type = line.type, point.size = point.size, point.shape = point.shape,
                       CI.color = CI.color, title = title, title.size= title.size, 
                       axis.label.size = axis.label.size, axis.tick.size = axis.tick.size, 
                       axis.x.label = axis.x.label,
                       axis.y.label = axis.y.label, 
                       facet.label.size = facet.label.size, facet.label.background = facet.label.background, scales = scales,
                       legend.title = legend.title,  legend.label = legend.label, legend.key.size = legend.key.size, legend.title.size = legend.title.size, 
                       legend.text.size = legend.text.size)
  }
  
}


#' @title Automatic Model Selection Results of IMU Object with Decomposed Processes
#' @description Creates a graph of the automatic model selection result with decomposed processes
#' @method autoplot auto.imu2
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
#' @param legend.title A \code{string} that indicates the title of legend. 
#' @param legend.label A \code{character vector} that indicates the labels on legend. 
#' @param legend.key.size A \code{double} that indicates the size of key (in centermeters) on legend. 
#' @param legend.title.size An \code{integer} that indicates the size of title on legend. 
#' @param legend.text.size An \code{integer} that indicates the size of key label on legend.
#' @param ... Additional options.
#' @author Wenchao
#' @return A panel containing the automatic model selection results of an IMU sensor.
autoplot.auto.imu2 = function(object, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                   line.type = NULL, point.size = NULL, point.shape = NULL,
                   CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                   axis.label.size = 13, axis.tick.size = 11, 
                   axis.x.label = expression(paste("Scale ", tau)),
                   axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                   facet.label.size = 13, facet.label.background = "#003C7D33", scales = "free_y", 
                   legend.title = '',  legend.label = NULL, legend.key.size = 1, legend.title.size = 13, 
                   legend.text.size = 13, ...){
  variable=variable2=value2=a_low=a_high=.x=NULL
  
  ###0. pre-process the object: auto.imu
  gmwm1 = object[[1]][[2]] #one gmwm object
  
  #what is num.sensor and ncols
  num.sensor = gmwm1$num.sensor
  ncols = sum(num.sensor)
  
  #what is freq
  freq = gmwm1$freq
  #what is alpha
  sig.level = gmwm1$alpha
  
  #what is axis
  axis = rep(0, ncols)
  sensor = rep(0, ncols)
  for(i in 1:ncols){
    axis[i] = object[[i]][[2]]$axis
    sensor[i] = object[[i]][[2]]$sensor
  }
  
  #Put data in the desired form
  obj.list = vector("list", ncols)
  
  #how large the data frame should be
  total.len = 0
  each.len = numeric(ncols)
  
  for (i in 1:ncols){
    obj.list[[i]] = object[[i]][[2]] 
    obj.list[[i]]$scales = obj.list[[i]]$scales/freq #freq conversion
    
    each.len[i] = length(obj.list[[i]]$wv.empir)
    total.len = total.len + each.len[i]
  }
  
  #initialize empty data frame with right number of rows
  #1) WV, low, high, theo 
  obj = data.frame(scales = numeric(total.len),
                   a_emp = numeric(total.len), 
                   a_low = numeric(total.len),
                   a_high = numeric(total.len),
                   z_theo = numeric(total.len),
                   axis = 'AXIS',
                   sensor = 'SENSOR', stringsAsFactors=FALSE)
  
  #2) decomped process
  all.model.desc = unlist( lapply(X = obj.list, FUN = function(x){addSpaceIfDuplicate(x$model$desc)} ) )
  # unique processes
  unique.model.desc = unique(all.model.desc)
  num.unique.model = length(unique.model.desc)
  
  num.col = num.unique.model + 3 # scales, process, axis, sensor
  process.df = as.data.frame(matrix(NA, ncol = num.col, nrow = total.len), stringsAsFactors = FALSE)
  colnames(process.df) = c('scales', unique.model.desc, 'axis', 'sensor')
  
  ###1. param checking
  if(CI){
    params = c('line.color', 'line.type', 'point.size', 'point.shape', 'CI.color', 'legend.label')
    requireLength = c(num.unique.model+3, num.unique.model+3, num.unique.model+3, num.unique.model+3, 1, num.unique.model+3)
    default = list(NULL, NULL,  NULL, NULL, "#003C7D", NULL)
    nullIsFine = c(rep(T,6))
  }else{
    params = c('line.color', 'line.type', 'point.size', 'point.shape', 'legend.label')
    requireLength = c(num.unique.model+2, num.unique.model+2, num.unique.model+2, num.unique.model+2, num.unique.model+2)
    default = list(NULL, NULL,  NULL, NULL, NULL)
    nullIsFine = c(rep(T,5))
  }
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  ###2. construct data frame
  t = 1
  for (i in 1:ncols){
    d = each.len[i]
    
    #1) WV, low, high, theo 
    obj[t:(t+d-1),] = data.frame(scales = obj.list[[i]]$scales, 
                                 a_emp = obj.list[[i]]$wv.empir,
                                 a_low = obj.list[[i]]$ci.low,
                                 a_high = obj.list[[i]]$ci.high,
                                 z_theo = obj.list[[i]]$theo,
                                 axis = axis[i], 
                                 sensor = sensor[i],
                                 stringsAsFactors=FALSE)
    
    #2) decomped process
    model.desc = c('scales', addSpaceIfDuplicate(obj.list[[i]]$model$desc), 'axis', 'sensor')
    process.df[t:(t+d-1), model.desc] = data.frame(obj.list[[i]]$scales, 
                                                   as.data.frame(obj.list[[i]]$decomp.theo),
                                                   axis = axis[i], 
                                                   sensor = sensor[i],
                                                   stringsAsFactors=FALSE)
    
    t = t +d
  }
  
  # 1) WV, low, high, theo
  if(CI){
    wv.df = data.frame(scales = obj$scales,
                     a_emp = obj$a_emp,
                     a_low = obj$a_low,
                     a_high = obj$a_high,
                     z_theo = obj$z_theo,
                     axis = obj$axis,
                     sensor = obj$sensor, stringsAsFactors = F)
    
    #construct the data frame to plot CI
    ci.df = data.frame(scales = obj$scales,
                        a_low = obj$a_low,
                        a_high = obj$a_high,
                        axis = obj$axis,
                        sensor = obj$sensor, stringsAsFactors = F)
    
  }else{
    wv.df = data.frame(scales = obj$scales,
                     a_emp = obj$a_emp,
                     z_theo = obj$z_theo,
                     axis = obj$axis,
                     sensor = obj$sensor, stringsAsFactors = F)
  }
  melt.wv.df = melt(wv.df, id.vars = c('scales', 'axis', 'sensor'))
  
  # 2) Process Decomp
  melt.process.df = melt(process.df, id.vars = c('scales', 'axis', 'sensor'), 
                         variable.name = 'variable2', value.name = 'value2')
  
  ###3. Auto-select the param
  # order to change the aethetics:
  # a_emp, a_high, a_low, decomposed prcesses, z_theo(sum)
  
  process.label = c(unique.model.desc, bquote("Implied WV "~nu*"("*hat(theta)*")"))
  process.color = getColors('Set1', num.unique.model, rm = 2) #remove blue color
  L = num.unique.model + 1 #decomposed process + theo
  
  if(CI == T){
    #CI.color
    if(is.null(CI.color)){CI.color = "#003C7D"}
    trans.CI.color = alpha(CI.color, transparence)
    
    #line type
    if(is.null(line.type)){
      line.type = c('solid','dotted', rep('solid', L) )}
    if(length(line.type) == (L+2) ){
      line.type = append(line.type, values = line.type[2], 2)}
    
    #line color
    if(is.null(line.color)){
      line.color = c("#003C7D", "#003C7D", process.color, "#F47F24")}
    if(length(line.color)== (L+2)){
      line.color = append(line.color, values = line.color[2], 2)}
    
    #point size
    if(is.null(point.size)){
      point.size = rep(0, L+2) }
    if(length(point.size)== (L+2)){
      point.size = append(point.size, values = point.size[2], 2)}
    
    #point shape
    if(is.null(point.shape)){
      point.shape = c(20, 20, rep(20,L-1), 1) }
    if(length(point.shape)== (L+2)){
      point.shape = append(point.shape, values = point.shape[2], 2)}
    
    #legend label
    if(is.null(legend.label)){
      legend.label = c(bquote("Empirical WV "~hat(nu)), 
                       bquote("CI("*hat(nu)*", "*.(1-sig.level)*")" ),
                       process.label) 
    }
    
    breaks = c('a_emp','a_low', unique.model.desc, 'z_theo')
    legend.fill = c(NA, trans.CI.color, rep(NA, L) )
    legend.linetype = c(line.type[1], 'blank', line.type[4:length(line.type)])
    legend.pointshape = c(point.shape[1], NA, point.shape[4:length(point.shape)])
    
  }else{
    #line.type
    if(is.null(line.type)){line.type = rep('solid', L+1) }
    
    #line color
    if(is.null(line.color)){line.color = c("#003C7D", process.color, "#F47F24")}
    
    #point size
    if(is.null(point.size)){point.size = rep(0, L+1)}
    
    #point shape
    if(is.null(point.shape)){point.shape =c(20, rep(20,L-1), 1) }
    
    #legend.label
    if(is.null(legend.label)){legend.label = c(expression(paste("Empirical WV ", hat(nu))), 
                                               process.label)}
    
    breaks = c('a_emp', unique.model.desc, 'z_theo')
  }
  
  ### 4. start to plot
  p = ggplot() +
    geom_line(data = melt.wv.df, mapping = aes(x = scales, y = value, linetype = variable, color = variable), na.rm = T) +
    geom_point(data = melt.wv.df, mapping = aes(x = scales, y = value, size = variable, shape = variable, color = variable), na.rm = T) +
    
    geom_line(data = melt.process.df, mapping = aes(x = scales, y = value2, linetype = variable2, color = variable2), na.rm = T) +
    geom_point(data = melt.process.df, mapping = aes(x = scales, y = value2, size = variable2, shape = variable2, color = variable2), na.rm = T) +
    
    scale_linetype_manual(name = legend.title, values = c(line.type), breaks = breaks, labels = legend.label ) +
    scale_shape_manual(name = legend.title, values = c(point.shape), breaks = breaks, labels = legend.label) +
    scale_size_manual(name = legend.title, values = c(point.size), breaks = breaks, labels = legend.label) +
    scale_color_manual(name = legend.title,values = c(line.color), breaks = breaks, labels = legend.label)
  
  if(CI){
    p = p + 
      geom_ribbon(data = ci.df, mapping = aes(x = scales, ymin = a_low, ymax = a_high), fill = trans.CI.color, show.legend = T, na.rm = T) + 
      guides(colour = guide_legend(override.aes = list(fill = legend.fill, linetype = legend.linetype, shape = legend.pointshape)))
  }
  
  if( background == 'white'){
    p = p + theme_bw() 
  }
  
  #what is y.lim
  y.lim.low = min(
    sapply(obj.list, FUN = function(x){
      0.8* min( c(x$ci.low, x$wv.empir) )
    })
  ) 
  y.lim.high = max(
    sapply(obj.list, FUN = function(x){
      1.05*max( c(x$wv.empir, x$ci.high))
    })
  )
  y.lim = c(y.lim.low, y.lim.high)
  
  p = p + facet_grid(sensor ~ axis, scales = scales) +
    
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    coord_cartesian(ylim = y.lim) +
    # set y lim makes scaels useless at all ==
    
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

#' @title Automatic Model Selection Results of IMU Object without Decomposed Processes
#' @description Creates a graph of the automatic model selection result without decomposed processes
#' @method autoplot auto.imu1
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
autoplot.auto.imu1 = function(object, CI = TRUE, background = 'white', transparence = 0.1, line.color = NULL, 
                             line.type = NULL, point.size = NULL, point.shape = NULL,
                             CI.color = "#003C7D", title = "Automatic Model Selection Results", title.size= 15, 
                             axis.label.size = 13, axis.tick.size = 11, 
                             axis.x.label = expression(paste("Scale ", tau)),
                             axis.y.label = expression(paste("Wavelet Variance ", nu)), 
                             facet.label.size = 13, facet.label.background = "#003C7D33",
                             scales = "free_y",...){
  value = variable = low = high = .x = NULL
  
  ###0. param checking
  
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
  
  checkParams(params = params, require.len = requireLength, default = default, null.is.fine = nullIsFine)
  
  ###1. pre-process the object: auto.imu
  
  #what is num.sensor and ncols
  num.sensor = object[[1]][[2]]$num.sensor
  ncols = sum(num.sensor)
  
  #what is freq
  freq = object[[1]][[2]]$freq
  
  #what is axis
  axis = rep(0, ncols)
  sensor = rep(0, ncols)
  for(i in 1:ncols){
    axis[i] = object[[i]][[2]]$axis
    sensor[i] = object[[i]][[2]]$sensor
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
  
  t = 1
  for (i in 1:ncols){
    d = each.len[i]
    
    obj[t:(t+d-1),] = data.frame(scales = obj.list[[i]]$scales/freq, # freq conversion
                                 emp = obj.list[[i]]$wv.empir,
                                 low = obj.list[[i]]$ci.low,
                                 high = obj.list[[i]]$ci.high,
                                 theo = obj.list[[i]]$theo,
                                 axis = axis[i], 
                                 sensor = sensor[i],
                                 stringsAsFactors=FALSE)
    t = t +d
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
      geom_ribbon(data = obj.CI, mapping = aes(x = scales, ymin = low, ymax = high), fill = alpha(CI.color, transparence), show.legend = F)
    
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
