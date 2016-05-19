context("Process WV Haar 1st Derivative - Unit Tests")

# Create a fake "zero" length vector
a = numeric()

# Overall taus
tau = 2^(1:4)

test_that("WN Derivative",{
  
  expect_equal(deriv_wn(tau), as.matrix( 1 / tau) )
  
})

test_that("RW Derivative",{
  
  expect_equal(deriv_rw(tau), as.matrix( (tau^2 + 2) / (12*tau) ) )
  
})

test_that("DR Derivative",{
  
  # Param value
  omega = .36
  
  expect_equal(deriv_dr(omega,tau), as.matrix((tau^2 * omega)/8))
  
})

test_that("QN Derivative",{
  
  expect_equal(deriv_qn(tau), as.matrix((6) / (tau^2)))
  
})

test_that("AR1 Derivative", {
  
  # Test Overalls
  phi = .23
  sigma2 = .2
  
  dphi = (2*sigma2)/((phi-1)^4*(phi+1)^2 * tau^2)*((phi^2-1)*tau*(-2*phi^(tau/2)+phi^(tau) - phi - 1) - (phi*(3*phi+2)+1)*(-4*phi^(tau/2)+phi^(tau)+3))
  dsigma2 = ((phi^2-1)*tau+2*phi*(-4*phi^(tau/2) + phi^(tau) + 3))/((phi-1)^3*(phi+1)*tau^2)
  
  testmat = cbind(dphi, dsigma2)
  
  expect_equal(deriv_ar1(phi, sigma2, tau), testmat, check.attributes = F)
})


test_that("MA1 Derivative", {
  
  # Test Overalls
  theta = .23
  sigma2 = .2
  
  dtheta = (sigma2*( 2*(theta+1)*tau - 6)) / (tau^2)
  dsigma2 = ((theta+1)^2*tau-6*theta)/(tau^2)
  
  testmat = cbind(dtheta, dsigma2)
  
  expect_equal(deriv_ma1(theta, sigma2, tau), testmat, check.attributes = F)
})

## Add in when ARMA(1,1) derivative implemented... 
# test_that("ARMA / ARMA(1,1) to WV", {
#   
#   # Test Overalls
#   phi = .62
#   theta = .23
#   sigma2 = .2
#
#   dphi    = ...
#   dtheta  = ...
#   dsigma2 = ...
#
#   testmat = cbind(dphi, dtheta, dsigma2)
#
#   # Test using deriv_arma11() to outside build
#   expect_equal(deriv_arma11(phi, theta, sigma2,tau), deriv_ma1(theta,sigma2,tau))
# 
#   # Test using deriv_ma1 or deriv_ar1 with deriv_arma11()
#   expect_equal(deriv_arma11(0, theta, sigma2,tau), deriv_ma1(theta,sigma2,tau))
#   expect_equal(deriv_arma11(phi, 0, sigma2,tau), deriv_ar1(phi,sigma2,tau))
#   
# })