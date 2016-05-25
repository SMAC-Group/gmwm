context("SARMA Object - Unit Tests")

# Create a fake "zero" length vector
a = numeric()

# Overall taus
tau = 2^(1:4)



ar  = 0.5753
i   = 0
ma  = -0.2709
sar = c(-0.2153, -0.2800)
si  = 0
sma = -0.4968
s   = 12
sigma2 = 1

# Model
mod = SARIMA(ar = ar, i = i, ma = ma, sar = sar, si = si, sma = sma, s = s,
             sigma2 = sigma2)

obj.desc = mod$obj.desc[[1]]

test_that("SARIMA object desc function",{
  
  expect_equal(sarma_objdesc(ar,ma,sar,sma,s,i,si), as.matrix( obj.desc ))
  
})

test_that("SARMA seasonal padding function",{
  
  
  # Number of seasonal components
  
  np  = obj.desc[1]
  nq  = obj.desc[2]
  nsp = obj.desc[3]
  nsq = obj.desc[4]
  ns  = obj.desc[6]
  
  p = np + ns * nsp;
  q = nq + ns * nsq;
  
  r_comp = matrix(c(p,q), ncol = 1)
  
  cpp_comp = sarma_calculate_spadding(np, nq, nsp, nsq, ns)
  
  expect_equal(cpp_comp, r_comp)
})

test_that("SARMA Components",{
  
  # Grab values
  np  = obj.desc[1]
  nq  = obj.desc[2]
  nsp = obj.desc[3]
  nsq = obj.desc[4]
  ns  = obj.desc[6]
  
  # Calculate total vec length
  p = np + ns * nsp;
  q = nq + ns * nsq;
  
  # Merge together
  r_comp = matrix(c(np,nq, nsp, nsq, ns, p,q), ncol = 1)
  
  # C++ 
  cpp_comp = sarma_components(obj.desc)
  
  # Test SARMA components C++ vs. R
  expect_equal(cpp_comp, r_comp)
  
})

test_that("SARMA Parameter Construction",{
  
  # Test Full (All vectors have value)
  r_comp = matrix(c(ar,ma,sar,sma),ncol = 1)

  cpp_comp = sarma_params_construct(ar,ma,sar,sma)
  
  expect_equal(r_comp, cpp_comp)
  
  # Test Empty (One vector does not have value)
  r_comp = matrix(c(ar,ma,sma),ncol = 1)
  
  cpp_comp = sarma_params_construct(ar,ma,numeric(),sma)
  
  expect_equal(r_comp, cpp_comp)
})

