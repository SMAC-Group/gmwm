#' @title Create Combination Matrix
#' @description Generates a matrix containing all the combination listings
#' @param n The number of variables.
#' @details Port expand.grid to C++ at a later time...
#' @return Returns a binary matrix (e.g. 1 or 0) entries
comb.mat = function(n){
  c = rep(list(1:0), n)
  expand.grid(c)
}

#' @title TS Model Checks
#' @description Stops the process in R if there is an issue with the model desc
#' @param desc A \code{character} vector containing \code{ts.model} terms from desc. 
#' @details Checks if there are two or more objects found of type: DR, QN, RW, or WN. In addition, it currently forbids ARMA Models.
#' @return Returns nothing if there is no issues. Otherwise, it will kill the process. 
select.desc.check = function(desc){
  models.active = count_models(desc)
  
  print(models.active)
  # Identifiability issues
  if(any( models.active[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, one of the supplied models will have identifiability issues. Please submit a new model.")
  }
  
  if(models.active[c("ARMA")] > 0){
    stop("Model selection with ARMA terms is NOT supported currently.")
  }
  
}

#' @title Automatically select appropriate model for a set of models
#' @description Runs through a model selection algorithm to determine the best model in a given set
#' @param ... Multiple \code{ts.model} objects
#' @param data A \code{vector}, \code{matrix}, \code{data.frame}, or \code{imu} object with either 1, 3, or 6 columns. 
#' @param bootstrap A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param alpha A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust A \code{boolean} that indicates whether to use robust estimation.
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B A \code{integer} that contains the amount of bootstrap replications
#' @param G A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @details The models MUST be nested within each other. If the models are not nested, the algorithm creates the "common denominator" model.
#' @return A \code{rank.models} object.
rank.models = function(..., data, bootstrap = F, alpha = 0.05, robust = F, eff = 0.6, B = 20, G = 1000){
  
  obj_list = list(...)
  numObj = length(obj_list)
  
  desc = vector("list", numObj) 
  
  for(i in 1:numObj){
    mod = obj_list[[i]]
    
    if(!is(mod, 'ts.model')){
      stop("Ill-formed ... request. Detected non-ts.model object in ... Please specify parameters with names")
    }
    select.desc.check(mod$desc)
    
    desc[[i]] = mod$desc
  }
  
  out = .Call('GMWM_find_full_model', PACKAGE = 'GMWM', x = desc)
  
  print(out)
  #out  = .Call('GMWM_auto_imu', PACKAGE = 'GMWM', )
  return(out)
}

#nested.models 

#' @title Automatically select appropriate model for IMU
#' @description Runs through a model selection algorithm to determine the best model
#' @param model A \code{ts.model} object that is the largest model to be tested.
#' @param data A \code{vector}, \code{matrix}, \code{data.frame}, or \code{imu} object with either 1, 3, or 6 columns. 
#' @param bootstrap A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param alpha A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust A \code{boolean} that indicates whether to use robust estimation.
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B A \code{integer} that contains the amount of bootstrap replications
#' @param G A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @return A \code{auto.imu} object.
auto.imu = function(data, model = 3*AR1()+WN()+RW()+QN()+DR(), bootstrap = F, alpha = 0.05, robust = F, eff = 0.6, B = 20, G = 1000){
  
  data = as.matrix(data)
  full.str = model$desc
  m = as.matrix(comb.mat(length(full.str)))
  m = m[-nrow(m),]
  
  out = .Call('GMWM_auto_imu', PACKAGE = 'GMWM', data, combs=m, full_model=full.str, alpha, compute_v = "fast", model_type = "sensor", K=1, H=B, G, robust, eff, bootstrap)
  model.names = build_model_set(m, full.str) 
  
  colnames(out[[1]][[1]][[1]]) = c("Obj Fun", "Optimism", "Criterion", "GoF P-Value")
  rownames(out[[1]][[1]][[1]]) = sapply(model.names, FUN = paste0, collapse=" ")
  
  m = out[[1]][[1]][[1]]
  
  out[[1]][[1]][[1]] = m[order(m[,3]), ]
  
  return(out)
}