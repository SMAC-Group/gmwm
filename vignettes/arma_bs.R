library(GMWM)

B = 50
N = 1000

mar = list(
  c(0.8, 0.1),c(0.8, 0.1),
  c(0.8897, -0.4858),c(0.8897, -0.4858),
  c(0.4609, -0.2185, -0.3671),c(0.4609, -0.2185, -0.3671)
)

mma = list(
  c(0.3),c(0.3),
  c(-0.2279, 0.2488),c(-0.2279, 0.2488),
  c(0.3012, 0.0759, 0.7996),c(0.3012, 0.0759, 0.7996)
)

out = vector(mode = "list", length = length(mma))

for(m in 1:length(mma)){
  
  ar = mar[[m]]
  ma = mma[[m]]
  va = 1
  
  p = length(ar)
  q = length(ma)
  S = p + q + 1
  
  cat(paste0("Working on model ARMA(",p,",",q,") with", if((m%%2)) "out contamination" else " contamination"))
  
  cat("\nUsing the following parameters\n")
  
  cat("AR: \n")
  print(ar)
  
  cat("MA: \n")
  print(ma)
  
  
  model = ARMA(ar,ma,va); 
  
  model2 = ARMA(p,q)
  
  timing = matrix(NA, B, 6)
  
  mod.mle = matrix(NA, B, S)
  mod.gmwm.guess = matrix(NA, B, S)
  
  mod.gmwm.guess.r = matrix(NA, B, S)
  
  mod.gmwm.exact = matrix(NA, B, S)
  mod.gmwm.exact.r = matrix(NA, B, S)
  
  mod.obj.guess = numeric(B)
  
  mod.obj.exact = numeric(B)
  
  global = vector("list", B)
  
  for(i in 1:B){
    print(paste0("Starting iteration:", i))
    timing[i, 1] = system.time({
      data = gen.gts(model, N)
      if((m %% 2)){
        data$data[round(runif(round(N*.1,0),1,N))] = rnorm(N*.1, 0, sd=10)
      }
    })[3]  
    
    timing[i, 2] = system.time({
      mle = arima(data$data,c(length(ar),0,length(ma)),include.mean=F)
    })[3]  
    
    print("Calculating MLE")
    mod.mle[i,] = c(mle$coef,mle$sigma2)
    print("Finished MLE")
    
    
    print("Calculating GMWM Guess")
    timing[i, 3] = system.time({
      est.guess = gmwm(model2, data, model.type="ssm", inference = F, compute.v='fast')
    })[3]
    
    mod.gmwm.guess[i,] = t(est.guess$estimate)
    print("Finished GMWM Guess")
    
    print("Calculating GMWM Guess Robust")
    timing[i, 4] = system.time({
      est.guess.r = gmwm(model2, data, model.type="ssm", inference = F, robust = T, compute.v='fast')
    })[3]
    
    mod.gmwm.guess.r[i,] = t(est.guess.r$estimate)
    
    
    print("Finished GMWM Guess Robust")
    
    print("Calculating Exact")
    timing[i, 5] = system.time({
      est.exact = gmwm(model, data, model.type="ssm", inference = F, compute.v='fast')
    })[3]
    
    mod.gmwm.exact[i,] = t(est.exact$estimate)
    print("Finished GMWM Exact")
    
    
    
    print("Calculating GMWM Exact Robust")
    timing[i, 6] = system.time({
      est.exact.r = gmwm(model, data, model.type="ssm", inference = F, robust = T, compute.v='fast')
    })[3]
    
    mod.gmwm.exact.r[i,] = t(est.exact.r$estimate)
    
    
    print("Finished GMWM Exact Robust")
    
    mod.obj.guess[i] = est.guess$obj.fun
    mod.obj.exact[i] = est.exact$obj.fun
    print(paste0("Finished Replication:  ", i))
    
  }

  obj = list(timing, 
             mod.mle, 
             mod.gmwm.guess,
             mod.gmwm.guess.r,
             mod.gmwm.exact,
             mod.gmwm.exact.r)
  out[[i]] = onj
}

p = 1
ar1 = cbind(mod.mle[,p],
            mod.gmwm.guess[,p],
            mod.gmwm.exact[,p],
            ifelse(mod.obj.exact > mod.obj.guess, mod.gmwm.guess[,p], mod.gmwm.exact[,p] ))



p = 2
ar2 = cbind(mod.mle[,p],
            mod.gmwm.guess[,p],
            mod.gmwm.exact[,p],
            ifelse(mod.obj.exact > mod.obj.guess, mod.gmwm.guess[,p], mod.gmwm.exact[,p] ))

p = 3
ma1 = cbind(mod.mle[,p],
            mod.gmwm.guess[,p],
            mod.gmwm.exact[,p],
            ifelse(mod.obj.exact > mod.obj.guess, mod.gmwm.guess[,p], mod.gmwm.exact[,p] ))



p = 4
sigma2 = cbind(mod.mle[,p],
               mod.gmwm.guess[,p],
               mod.gmwm.exact[,p],
               ifelse(mod.obj.exact > mod.obj.guess, mod.gmwm.guess[,p], mod.gmwm.exact[,p] ))


par(mfrow = c(2,3))
boxplot(ar1, col = "lightgrey")
abline(h = ar[1], lwd = 2, col = 2)

boxplot(ar2, col = "lightgrey")
abline(h = ar[2], lwd = 2, col = 2)

boxplot(ma1, col = "lightgrey")
abline(h = ma[1], lwd = 2, col = 2)

boxplot(sigma2, col = "lightgrey")
abline(h = va, lwd = 2, col = 2)

save(ar1,ar2,ma1,sigma2, file="NLM.rda")


system.time({
  est.param = gmwm(model, data)
})



summary(est.param)

summary(est.guess)

summary(est.guess2)

mm = list(ar = ar, ma = ma )
