context("Process to WV Haar - Unit Tests")

# Create a fake "zero" length vector
a = numeric()


# Overall taus
tau = 2^(1:4)

test_that("WN to WV",{
  
  # Param value
  sigma2 = .01
  
  expect_equal(wn_to_wv(sigma2,tau), as.matrix( sigma2 / tau))
  
})

test_that("RW to WV",{
  # Param value
  gamma2 = .01
  
  expect_equal(rw_to_wv(gamma2,tau), as.matrix( (gamma2 * (tau^2 + 2.0)) / (12.0*tau)))
  
})

test_that("DR to WV",{
  
  # Param value
  omega = .36
  
  expect_equal(dr_to_wv(omega,tau), as.matrix(tau^2 * omega^2 / 16))
  
})

test_that("QN to WV",{
  
  # Param value
  q2 = .36
  
  expect_equal(qn_to_wv(q2,tau), as.matrix(6*q2/tau^2))
  
})

test_that("AR1 to WV", {
  
  # Test Overalls
  phi = .23
  sigma2 = .2

  # Write a quick R test
  ar1_to_wv_verR = function(phi,sigma2, tau){
    as.matrix(
      (sigma2*(2*phi*(3 - 4*phi^(tau/2) + 
                      phi^tau) + (-1 + phi^2)*
               tau))/((-1 + phi)^3*(1 + phi)*
                        tau^2)
    )
  }
  
  # Test C++ with R version
  expect_equal(ar1_to_wv(phi,sigma2,tau),ar1_to_wv_verR(phi,sigma2,tau))
  
  # Test C++ with general ARMA(1,0)
  expect_equal(ar1_to_wv(phi,sigma2,tau), arma_to_wv(phi,a,sigma2,tau) )
})


test_that("MA1 to WV", {
  
  # Test Overalls
  theta = .23
  sigma2 = .2
  
  expect_equal(ma1_to_wv(theta,sigma2,tau), arma_to_wv(a,theta,sigma2,tau) )
})

test_that("ARMA / ARMA(1,1) to WV", {
  
  # Test Overalls
  phi = .62
  theta = .23
  sigma2 = 1
  
  # Quick R implementation
  arma11_to_wv_r = function(phi, theta, sigma2, tau){
    as.matrix(
      (-2*sigma2*((-(theta + phi))*(1 + theta*phi)*
                    (3 - 4*phi^(tau/2) + phi^tau) - 0.5*(1 + theta)^2*(-1 + phi^2)*tau)) /
        ((-1 + phi)^3*(1 + phi)*tau^2)
    )
  }
  

  # Test C++ vs. R
  expect_equal(arma11_to_wv(phi,theta,sigma2,tau), arma11_to_wv_r(phi,theta,sigma2,tau) )
  
  # Test C++ vs. Generic C++
  expect_equal(arma11_to_wv(phi,theta,sigma2,tau), arma_to_wv(phi,theta,sigma2,tau) )
  
  # Test individual processes
  expect_equal(arma11_to_wv(0, theta, sigma2,tau), ma1_to_wv(theta,sigma2,tau))
  expect_equal(arma11_to_wv(phi, 0, sigma2,tau), ar1_to_wv(phi,sigma2,tau))

})


test_that("Process to Theoretical Sum WV", {
  
  # Settings
  phi = .62
  phi_sigma2 = 1
  wn = 2
  dr = 0.0001
  qn = 0.2
  rw = 0.025
  theta = 0.1
  theta_sigma2 = 0.5
  
  phi2 = 0.78
  theta2 = 0.32
  sigma2 = 1.25
  
  phi2arma = 0.21
  theta2arma = 0.84
  sigma2arma = 1.5
  
  # Separate models
  wv.theo.split.r = cbind(ar1_to_wv(phi,phi_sigma2, tau),
                          ma1_to_wv(theta,theta_sigma2,tau),
                          arma11_to_wv(phi2,theta2,sigma2,tau),
                          arma_to_wv(phi2arma,theta2arma,sigma2arma,tau),
                          wn_to_wv(wn, tau),
                          dr_to_wv(dr, tau),
                          qn_to_wv(qn, tau),
                          rw_to_wv(rw, tau)
                          )
  # Total WV
  wv.total.r = as.matrix(rowSums(wv.theo.split.r))
  
  # Test Overalls
  model = AR1(phi, phi_sigma2) +  MA1(theta, theta_sigma2) + 
              ARMA11(phi2, theta2, sigma2) + 
              ARMA(phi2arma,theta2arma,sigma2arma) +
              WN(wn) + DR(dr) + QN(qn) + RW(rw) 
  
  
  wv.theo.split = decomp_theoretical_wv(model$theta, model$desc, model$obj.desc, tau)
  
  # Test Split Decomposition of WV C++ vs. R
  expect_equal(wv.theo.split, wv.theo.split.r)
  
  # Test Row Sum Split Decomposition vs. R
  expect_equal(decomp_to_theo_wv(wv.theo.split),  wv.total.r)
  
  # Test Individual Process Summation C++ vs. R
  expect_equal(theoretical_wv(model$theta, model$desc, model$obj.desc, tau), wv.total.r)
  
})