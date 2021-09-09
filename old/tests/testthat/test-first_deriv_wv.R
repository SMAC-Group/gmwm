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
  
  dphi = (2*sigma2*((-(3 - 4*phi^(tau/2) + phi^tau))*(1 + phi*(2 + 3*phi)) + 
                    (-1 + phi^2)*(-1 - phi - 2*phi^(tau/2) + phi^tau)*tau)) /
         ((-1 + phi)^4*(1 + phi)^2*tau^2)
  
  dsigma2 = ((phi^2 - 1)*tau + 
              2*phi*(phi^tau - 4*phi^(tau/2) + 3)) /
            ((phi - 1)^3*(phi + 1)*tau^2)
  
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

test_that("ARMA / ARMA(1,1) to WV", {

  # Test Overalls
  phi    = 0.23
  theta  = 0.31
  sigma2 = 0.20

  # R Coded Values
  dphi    = (1/((-1 + phi)^4*(1 + phi)^2*tau^2))*2*sigma2*((-(3 - 4*phi^(tau/2) + phi^tau))*(1 + phi*(2 + 3*phi) + theta^2*(1 + phi*(2 + 3*phi)) + 2*theta*(1 + phi*(3 + phi + phi^2))) + 
              ((-(1 + theta)^2)*(-1 + phi)*(1 + phi)^2 - 
                 2*phi^(tau/2 - 1)*(theta + phi)*
                 (1 + theta*phi)*(-1 + phi^2) + 
                 phi^(tau - 1)*(theta + phi)*
                 (1 + theta*phi)*(-1 + phi^2))*tau)
  dtheta  = (2*sigma2*((1 + 2*theta*phi + phi^2)*(3 - 4*phi^(tau/2) + phi^tau) +
                         (1 + theta)*(-1 + phi^2)*tau)) / 
                       ((-1 + phi)^3*(1 + phi)*tau^2)
  dsigma2 = ((-2*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau/2) + phi^tau) - 
                   (1/2)*(1 + theta)^2*(-1 + phi^2)*tau))/
              ((-1 + phi)^3*(1 + phi)*tau^2))

  testmat = cbind(dphi, dtheta, dsigma2)

  # Test deriv_arma11 using deriv_ar1
  expect_equal(deriv_arma11(phi, 0, sigma2,tau)[,-2], deriv_ar1(phi,sigma2,tau))
  
  # Test deriv_arma11 using deriv_ma
  expect_equal(deriv_arma11(0, theta, sigma2, tau)[,-1], deriv_ma1(theta,sigma2,tau))
  
  # Test deriv_arma11 with R implementation
  expect_equal(deriv_arma11(phi, theta, sigma2, tau), testmat, check.attributes = F)
  

})
