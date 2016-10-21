# Copyright (C) 2014 - 2016  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Place Legend
#' @description This function decides where the legend should be put (top left or bottom left)
#' @param wv_1 A \code{double} that indicates the first value of \code{wv.empir}
#' @param low_n A \code{doble} that indicates the last value of \code{ci.low}
#' @param high_n A \code{dobule} that indicates the last value of \code{ci.high}
#' @return A numeric vector containing 4 elements. The first two elements indicate legend.justification, the last two elements indicate legend.position (see \code{?theme}).
#' @keywords internal
placeLegend = function(wv_1, low_n, high_n){
  if(log10(wv_1) > ( log10(low_n) + log10(high_n) )/2 ){
    # legend should be placed in bottom left
    legend.justification = c(0,0)
    legend.position = c(0,0)
    #x = wv_1[1]/xlim_length
    #y = high_n/ylim_length
  }
  else{
    # legend should be placed in top left
    legend.justification = c(0,1)
    legend.position = c(0,1)
    #x = wv_1[1]/xlim_length
    #y = low_n/ylim_length
  }
  return( c(legend.justification, legend.position) )
  
}

#' @title Frequent Graph Setting for Paper
#' @description This function sets some parameters such as plot.margin.
#' @return A ggplot2 panel containing the frequent graph setting for paper.
#' @keywords internal
paperSetting = function(){
  
  p = theme(axis.title.y=element_text(margin=margin(0,22.5,0,0)), 
            axis.title.x=element_text(margin=margin(22.5,0,0,0)), 
            plot.title = element_text(margin=margin(0,0,20,0)), 
            plot.margin=unit(c(0.7,0.1,0.7,0.7),"cm"))
  
  return(p)
}

#' @title Emulate ggplot2 default color palette
#' @description Autogenerate a colors according to the ggplot selection mechanism. 
#' @param n An \code{integer} indicating how many colors user wants.
#' @return A \code{vector} containing \code{n} colors
#' @author John Colby
#' @keywords internal
ggColor <- function(n) {
  hues = seq(15, 375, length=n+1)
  rev(hcl(h=hues, l=70, c=100)[1:n])
}

#' @title Get the model in a \code{gmwm} object
#' @description Extracts and formats the model string.
#' @param object A \code{gmwm} object
#' @return A \code{string} containing the model
#' @keywords internal
getModel.gmwm = function(object){
  if( !is(object, 'gmwm') ){
    stop('It must be a gmwm object')
  }
  model.desc = object$model.hat$desc
  count.map = count_models(model.desc)
  all.model = names(count.map)
  
  model = ''
  for (i in 1:length(all.model)){
    if(length(model)==1){
      if(count.map[all.model[i]] == 1){
        model = bquote(.(all.model[i]))
      }else if(count.map[all.model[i]]>1){
        model = bquote(.(count.map[all.model[i]])%*%.(all.model[i]))
      }
    }else{
      if(count.map[all.model[i]] == 1){
        model = bquote(.(model) * "+" *.(all.model[i]))
      }else if(count.map[all.model[i]]>1){
        model = bquote(.(model) * "+" *.(count.map[all.model[i]])%*%.(all.model[i]) )
      }
    }
  }
  return(as.expression(model) )
}


#' @title Order the Model
#' @description Orders the model and changes it to the correct format
#' @param models A vector of \code{string} that specifies the models
#' @details If the \code{models} are c("AR1", "WN", "AR1", "WN", "AR1+WN+AR1+WN"), it will be converted to 
#' c("AR1-1", "WN-1", "AR1-2", "WN-2", "AR1+WN+AR1+WN").
#' 
#' This function is used in \code{gen.lts()}
#' @keywords internal
#' @examples 
#' models = c("AR1", "WN", "AR1", "WN", "AR1+WN+AR1+WN")
#' new.models = orderModel(models)
#' new.models
#' 
#' models = c('AR1', 'QN', 'WN', 'AR1+QN+WN')
#' new.models = orderModel(models)
#' new.models
orderModel = function(models){
  count = table(models)
  if( any(count>1)){
    multi.models = names( count[count>1] )
    
    for(model in multi.models){
      num = count[model]
      models[models == model] = paste( rep(model, num), rep('-', num), 1:num, sep = '' )
    }
    
    return(models)
  }else{
    return(models)
  }
}

#' @title Get the Computation Method and Efficiency of \code{gmwm} Object
#' @description The computation method (classical/robust) and efficiency will be returned in a certain format.
#' @param x A \code{gmwm} object
#' @details Used in \code{compare.eff()}.
#' @keywords internal
#' @examples 
#' n = 1000
#' x = gen_gts(n, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' GMWM1 = gmwm(2*AR1()+RW(), data = x, robust  = TRUE, eff = 0.9)
#' GMWM2 = gmwm(2*AR1()+RW(), data = x, robust  = FALSE)
#' 
#' formatRobustEff(GMWM1)
#' formatRobustEff(GMWM2)
formatRobustEff = function(x){
  if(!x$robust){
    return('Classical')
  }else{
    res = paste0('Robust eff. ', x$eff)
    return(res)
  }
} 

#' @title Validity of the Object List for \code{compare.eff}
#' @description Check whether the object list is valid for function \code{compare.eff}.
#' @param obj.list A \code{list} of \code{gmwm} object
#' @details 
#' A valid object list should contain \code{gmwm} objects constructed by the same data
#' and the same model.
#' 
#' @keywords internal
is.validCompEffObj = function(obj.list){
  
  # All gmwm object
  sapply(obj.list, FUN = 
                     function(x){ 
                       if( !is.gmwm(x) ){
                          stop('The function can only work on gmwm object.')}
                     })
  
  #constructed by same data
  expectDiff = sapply(obj.list, FUN = 
                      function(x){x$expect.diff})
  if( any(expectDiff!=expectDiff[1]) ){
    stop('This function can only operate on models constrcuted by the same data.') 
  }
  
  #same model
  count.map = sapply(obj.list, FUN = function(x){
    count_models( x$model.hat$desc )
  })
  sameModel = apply(count.map, 2, identical, count.map[,1])
  
  if(any(sameModel==F)){
    stop('gmwm objects are not constructed by the same model.')
  }
  
}

#' @title Add Space to Avoid Duplicate Elements
#' @description Add space to every element if there are duplicates in the vector.
#' @param x A \code{character vector}
#' @keywords internal
#' @examples 
#' ##no duplicate
#' x1 = c('GMWM2', 'GMWM1', 'GMWM3')
#' addSpaceIfDuplicate(x1)
#' 
#' ##duplicate
#' x2 = c('GMWM3', 'GMWM4', 'GMWM3', 'GMWM4', 'GMWM5', 'GMWM6','GMWM3')
#' addSpaceIfDuplicate(x2)
addSpaceIfDuplicate = function(x){
  res = x
  count.table = table(x)
  
  if ( any( count.table > 1)  ){
    
    dup.item = count.table[ count.table>1 ]
    
    for(each in names(dup.item) ){
      index = which(x == each)
      
      for(i in 2:length(index)){
        res[ index[i]  ] = paste0(x[ index[i] ], paste0(rep(' ',times = (i-1) ), collapse = ''))
      }
     
    }
    
  }
  return(res)
}


#' @title Get \code{gmwm} Efficiency Values
#' @description Get efficiency values from a list of \code{gmwm} object
#' @param obj.list A \code{list} of \code{gmwm} object
#' @return A \code{numeric vector}.
#' @keywords internal
#' @details 
#' 
#' If the object is computed by classical method, it will return 1.1. The reason is:
#' 
#' It's possible for user to create one object by robust method with eff=1, though it is exactly same
#' as classical method. In this case, if we want the classical method to always appear on the 
#' top left corner,  number larger than 1 (e.g. 1.1) is used. This setting makes it easy to draw the graph.
#' 
#' @examples 
#' set.seed(8836)
#' n = 1000
#' x = gen_gts(n, AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1))
#' GMWM1 = gmwm(2*AR1()+RW(), data = x, robust  = FALSE)
#' GMWM2 = gmwm(2*AR1()+RW(), data = x, robust  = TRUE, eff = 0.1)
#' GMWM3 = gmwm(2*AR1()+RW(), data = x, robust  = TRUE, eff = 0.6)
#' obj.list = list(GMWM1, GMWM2, GMWM3)
#' getEff(obj.list)
getEff = function(obj.list){
  res = sapply(obj.list, FUN = function(x){
    
    if(x$robust){return(x$eff)
    }else{return(1.1)} #classical method is robust with eff=1
    #eff = Var(Robust)/Var(Classical)
  })
  
  return(res)
}

#' @title Get N Colors
#' @description Creates n colors from specific palette
#' @param palette A \code{string} that indicates the name of the palette.
#' @param n An \code{integer} that specifies the number of colors to be generated.
#' @param rm An \code{integer vector} that specifies the index of colors to be removed.
#' @keywords internal
#' @details 
#' Currently, only the palette 'Set1' is implemented in this function.
#' 
#' 'Set1' palette in the package \code{RColorBrewer} originally has 9 colors, but we delete the 
#' yellow color (the 6th one) since it is too much vibrant. 
#' 
#' If \code{n} is larger than the mamixum number of colors in the palette, the function
#' will repeat the palette until \code{n} colors is created.
#' 
#' @examples 
#' getColors('Set1', 10)
getColors = function(palette, n, rm = NULL){
  
  if(palette == 'Set1'){
    color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628" ,"#F781BF", "#999999")
  }
  
  if(!is.null(rm)){
    color = color[-rm]
  }
  
  modulus = n %/% length(color)
  remainder = n %% length(color)
  result = c( rep(color, times = modulus), color[1:remainder] )
  return(result)
}

#' @title Check the Parameters
#' @description Check the user supplied parameters and assign them to the default if they are wrong.
#' @param params A \code{character vector} that specifies the user supplied parameters.
#' @param require.len An \code{integer vector} that specifies the required length of each parameter.
#' @param default A \code{list} that specifies the default of each parameter.
#' @param null.is.fine A \code{boolean vector} to indicate whether \code{NULL} is fine for parameters.
#' @param env An \code{environment} to use.
#' @keywords internal
#' @details 
#' 
#' The user supplied parameters are usually \code{line.color}, \code{line.type}, \code{point.size}, 
#' \code{point.shape}, \code{CI.color} and \code{legend.label}.
#' 
#' This function will check whether the required length of the parameter is met. If not, it will assign the 
#' default value to that parameter.
#' 
checkParams = function(params, require.len, default, null.is.fine, env = parent.frame()){
  
  for (i in 1:length(params)){
    
    one_param = params[i]
    value = get(one_param, envir = env)
    
    if( length(value)!=require.len[i]){
      
      isNull = is.null(value)
      
      if( (!isNull) || (!null.is.fine[i]) ){
        
        warning(paste('Parameter', one_param, 'requires', require.len[i],'elements,','but', length(value),
                      'is supplied.','Default setting is used.'))
      }
      
      assign(one_param, default[[i]], envir = env)
    }
  }
  
}


#' @title Add Units to Sensor
#' @description Add corresponding units for gyroscope and accelerometer sensors.
#' @param units \code{NULL} or a two-element vector indicating the units of gyroscope and accelerometer sensors.
#' @param sensor A \code{character vector} including all sensors in imu objects.
#' @keywords internal
#' @details 
#' This function is used in \code{plot.wvar.imu} and \code{plot.auto.imu}.
addUnits = function(units, sensor){
  #add units to sensors
  if(!is.null(units)){
    
    sensor[sensor == "Gyroscope"] = as.character(as.expression( bquote('Gyroscope ' * .(units[[1]])) ))
    sensor[sensor == "Accelerometer"] = as.character(as.expression( bquote("Accelerometer " * .(units[[2]])) ))
  }
  return(sensor)
}

#' @title Check Equal Attributes in a GMWM List
#' @description Check whether the \code{gmwm} objects in \code{obj.list} have equal attributes.
#' @return A list of square matrix indicating whether the objects have the same attributes.
#' @param obj.list A list of \code{gmwm} objects.
#' @keywords internal
#' @author Wenchao
#' @details 
#' This function is used in \code{\link{compare.models}}.
checkEquality = function(obj.list){
  res = list(scales = NULL, ci.low = NULL, ci.high = NULL, wv.empir = NULL, theo = NULL)
  elements = names(res)
  
  numObj = length(obj.list)
  store = matrix(data = NA, nrow = numObj, ncol = numObj)
  
  for(element in elements){
    
    for(i in 1:numObj){
      for(j in 1:numObj){
        store[i,j] = identical(obj.list[[i]][[element]], obj.list[[j]][[element]])
      }
    }
    res[[element]] = store
  }
  
  res$bounds = res$ci.low & res$ci.high
  res$ci.low = NULL
  res$ci.high = NULL
  
  #If scales are different, everything is different
  if(any(res$scales == F)){
    res = lapply(X = res, FUN = function(x){ x & res$scales })
  }
  
  return(res)
}


#' @title Autofill A Vector
#' @description Autofill a vector to a target length with the specified element
#' @return A new \code{vector}
#' @param v A \code{vector}
#' @param len An \code{integer} indicating the target length
#' @param fillwith The specified element used to fill the empty cells
#' @keywords internal
#' @details 
#' If the length of \code{v} is less than \code{len}, the function will introduce 
#' \code{fillwith} for empty cells.
#' 
#' Nothing will be done if the length of \code{v} is equal to or larger than \code{len}.
#' 
autofill = function(v, len, fillwith = NA){
  v.len = length(v)
  
  if(v.len < len){
    v = c(v, rep(fillwith, len-v.len))
  }
  return(v)
}
