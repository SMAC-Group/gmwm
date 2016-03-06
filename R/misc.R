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
#' x = gen.gts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
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
#' x = gen.gts(AR1(phi = .1, sigma2 = 1) + AR1(phi = 0.95, sigma2 = .1), n)
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

