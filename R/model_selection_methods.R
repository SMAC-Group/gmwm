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
  
  # Identifiability issues
  if(any( models.active[c("DR","QN","RW","WN")] >1)){
    stop("Two instances of either: DR, QN, RW, or WN has been detected. As a result, one of the supplied models will have identifiability issues. Please submit a new model.")
  }
  
  if(models.active[c("ARMA")] > 0){
    stop("Model selection with ARMA terms is NOT supported currently.")
  }
  
}

#' @title Formats the model score matrix
#' @description The model score matrix receives the appropriate model numbering, col descriptors, and ordering of models.
#' @param out A \code{list} containing the model matrix
#' @param num.models A \code{integer} that indicates the number of models fit. 
#' @return An updated matrix in the ith position of the list. 
pretty.model.score = function(out, model.names){

  colnames(out) = c("Obj Fun", "Optimism", "Criterion", "GoF P-Value")
  rownames(out) = sapply(model.names, FUN = paste0, collapse=" ")
  
  #out = out[order(out[,3]), ]

  out
}


output.format = function(out, model.names, scales, N, alpha, robust, eff, B, G, seed){
  desc = model.names[out[[1]][[2]]]
  
  model.ts = desc.to.ts.model(desc[[1]])
  
  out[[1]] = pretty.model.score(out[[1]][[1]], desc)
  
  gmwm.obj = out[[2]]
  estimate = gmwm.obj[[1]]
  rownames(estimate) = model.ts$process.desc
  init.guess = gmwm.obj[[2]]
  rownames(init.guess) = model.ts$process.desc
  
  model.hat = model.ts
  
  model.hat$starting = F  
  
  model.hat$theta = as.numeric(estimate)
  
  # Release model
  out[[2]] = structure(list(estimate = estimate,
                                 init.guess = init.guess,
                                 wv.empir = gmwm.obj[[3]], 
                                 ci.low = gmwm.obj[[4]], 
                                 ci.high = gmwm.obj[[5]],
                                 orgV = gmwm.obj[[7]],
                                 V = gmwm.obj[[6]],
                                 omega = gmwm.obj[[12]],
                                 obj.fun = gmwm.obj[[11]],
                                 theo = gmwm.obj[[9]],
                                 decomp.theo = gmwm.obj[[10]],
                                 scales = scales, 
                                 robust = robust,
                                 eff = eff,
                                 model.type = "imu",
                                 compute.v = "fast",
                                 augmented = F,
                                 alpha = alpha,
                                 expect.diff = gmwm.obj[[8]],
                                 N = N,
                                 G = G,
                                 H = B,
                                 K = 1,
                                 model = model.ts, # ADD THIS IN AT A LATER TIME!
                                 model.hat = model.hat,
                                 starting = TRUE,
                                 seed = seed), class = "gmwm")
  
  out
}


#' @title Automatically select appropriate model for a set of models
#' @description Runs through a model selection algorithm to determine the best model in a given set
#' @param ... Multiple \code{ts.model} objects
#' @param data A \code{vector}, \code{data.frame}, \code{matrix}, or \code{gts} object with 1 column.
#' @param bootstrap A \code{bool} that is either true or false to indicate whether we use bootstrap or asymptotic By default, we use asymptotic.
#' @param alpha A \code{double} that indicates the level of confidence for the WV CI.
#' @param robust A \code{boolean} that indicates whether to use robust estimation.
#' @param eff A \code{double} between 0 and 1 that indicates the efficiency for the robust estimation.
#' @param B A \code{integer} that contains the amount of bootstrap replications
#' @param G A \code{integer} that indicates the amount of guesses for caliberating the startup.
#' @details The models MUST be nested within each other. If the models are not nested, the algorithm creates the "common denominator" model.
#' @return A \code{rank.models} object.
rank.models = function(data, models=list(AR1()+WN(), AR1()), nested = F, bootstrap = F, model.type="ssm", alpha = 0.05, robust = F, eff = 0.6, B = 50, G = 100000, seed = 1337){
  
  set.seed(seed)
  numObj = length(models)
  desc = vector("list", numObj) 
  
  for(i in 1:numObj){
    mod = models[[i]]
    
    if(!is(mod, 'ts.model')){
      stop("Ill-formed ... request. Detected non-ts.model object in ... Please specify parameters with names")
    }
    select.desc.check(mod$desc)
    
    desc[[i]] = mod$desc
  }
  
  
  if(nested == F){
    full.str = .Call('GMWM_find_full_model', PACKAGE = 'GMWM', x = desc)
    
    if(!any(sapply(desc, function(x, want) isTRUE(all.equal(x, want)),  full.str)) ){
      print("Creating a Common Denominator Model!")
      desc[[length(desc)+1]] = full.str
    }
    desc = vector_to_set(desc)
  }else{
    
    full.str = models[[1]]$desc
    
    m = as.matrix(comb.mat(length(full.str)))
    m = m[-nrow(m),]
    
    desc = build_model_set(m, full.str)
    
  }
  
  out = .Call('GMWM_rank_models', PACKAGE = 'GMWM', data, model_str=desc, full_model=full.str, alpha, compute_v = "fast", model_type = model.type, K=1, H=B, G, robust, eff, bootstrap)
  
  N = length(data)
  nlevels =  floor(log2(N))
  scales = .Call('GMWM_scales_cpp', PACKAGE = 'GMWM', nlevels)
  
  out[[1]] = output.format(out[[1]], desc, scales, N, alpha, robust, eff, B, G, seed)  
  
  class(out) = c("rank.models")
  
  out
}

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
auto.imu = function(data, model = 3*AR1()+WN()+RW()+QN()+DR(), bootstrap = F, alpha = 0.05, robust = F, eff = 0.6, B = 50, G = 100000, seed = 1337){
  
  #check object
  if(is.null(data) || !is(data,"imu") ) {
    stop('Object must an imu object via imu()')
  }
  
  # Extract data for IMU Injection
  data.in = as.matrix(data$data) 
  sensors = data$sensor
  num.sensor = data$num.sensor
  axis = data$axis

  # Set seed for reproducible results
  # Need to figure out a way to do set same seed on each model generation.
  set.seed(seed)
  
  # Set up data and models for automatic processing
  full.str = model$desc
  m = as.matrix(comb.mat(length(full.str)))
  m = m[-nrow(m),]
  
  out = .Call('GMWM_auto_imu', PACKAGE = 'GMWM', data.in, combs=m, full_model=full.str, alpha, compute_v = "fast", model_type = "sensor", K=1, H=B, G, robust, eff, bootstrap)
  
  # Handle post processing
  
  # Get model names
  model.names = build_model_set(m, full.str) 
  
  # Get basic data info
  N = nrow(data.in)
  nlevels =  floor(log2(N))
  scales = .Call('GMWM_scales_cpp', PACKAGE = 'GMWM', nlevels)
  
  # Get set statements for Wenchao
  n.gyro = num.sensor[1]
  n.acc = num.sensor[2]
  
  a.gyro = 0
  a.acc = 0
  
  for(i in 1:ncol(data.in)){
    obj = out[[i]]
    obj = output.format(obj, model.names, scales, N, alpha, robust, eff, B, G, seed)
    
    obj.gmwm = obj[[2]]
    if(n.acc != 0 && a.acc != (n.acc - 1)){
      a.acc = a.acc + 1 
      obj.gmwm$sensor = "Gyroscope"
      obj.gmwm$axis = axis[a.acc]
    } else if(n.gyro != 0 && a.gyro != (n.gyro - 1)){
      a.gyro = a.gyro + 1 
      obj.gmwm$sensor = "Accelerometer"
      obj.gmwm$axis = axis[a.gyro]
    }
    
    obj.gmwm$num.sensor = num.sensor
    
    obj[[2]] = obj.gmwm
    out[[i]] = obj
  }
  
  class(out) = c("auto.imu","rank.models")
  out  
}

print.rank.models = function(object, ...){
  summary.rank.models(object)
}

print.auto.imu = function(object, ...){
  summary.auto.imu(object)
}

summary.auto.imu = function(object, digits = 4, ...){
  
  n.process = length(object)
  
  cat(paste0("There were ", n.process, " observed\n\n"))
  
  for(i in 1:n.process){
    out = object[[i]][[1]]
    cat(paste0("The model ranking for data column ", i, ": \n"))
    rownames(out) = paste0(1:nrow(out), ". ", rownames(out) )
    
    print(round(out,digits))
    
    cat("\n")
  }
}

summary.rank.models = function(object, digits = 4, ...){
  cat("The model ranking is given as: \n")
  out = object[[1]][[1]]
  
  rownames(out) = paste0(1:nrow(out), ". ", rownames(out) )
  
  print(round(out,digits))
}