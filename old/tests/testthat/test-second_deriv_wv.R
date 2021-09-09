context("Process WV Haar 2nd Derivative - Unit Tests")

# Create a fake "zero" length vector
a = numeric()

# Overall taus
tau = 2^(1:4)

test_that("DR 2nd Derivative",{
  
  # Param value
  omega = .36
  
  expect_equal(deriv_2nd_dr(tau), as.matrix((tau^2/8)))
  
})

test_that("AR1 2nd Derivative", {
  # Test Overalls
  phi    = .23 # 0  
  sigma2 = .2
  
  d2phi = (1/((phi - 1)^5*(phi + 1)^3* tau^2))*
    (2*sigma2*((phi^2 - 1)*tau*(2*((7*phi + 4)*phi + 1)*phi^(tau/2 - 1) - 
                                  ((7*phi + 4)*phi + 1)*phi^(tau - 1) + 
                                  3*(phi + 1)^2) + 
                                  (phi^2 - 1)^2*tau^2*(phi^(tau/2) - 1)*phi^(tau/2 - 1) + 
                 4*(phi^2 + phi + 1)*(3*phi + 1)*(phi^tau - 4*phi^(tau/2) + 3)))
  
  dphisigma = (2*((-(3 - 4*phi^(tau/2) + phi^tau))*(1 + phi*(2 + 3*phi)) + (-1 + phi^2)*
                    (-1 - phi - 2*phi^(tau/2) + phi^tau)*tau))/((-1 + phi)^4*(1 + phi)^2*tau^2)
  
  dsigma4 = 0

  testmat = cbind(d2phi, dphisigma,dsigma4)
  
  expect_equal(deriv_2nd_ar1(phi, sigma2, tau), testmat, check.attributes = F)
})


test_that("MA1 2nd Derivative", {
  
  # Test Overalls
  theta = .23
  sigma2 = .2
  
  d2theta = 2 * sigma2 / tau
  dthetasigma2 = (-6 + 2*(1 + theta)*tau)/tau^2
  d2sigma2 = 0
  
  testmat = cbind(d2theta, dthetasigma2, d2sigma2)
  
  expect_equal(deriv_2nd_ma1(theta, sigma2, tau), testmat, check.attributes = F)
})


test_that("ARMA11 2nd Derivative", {

  # Test Overalls

  phi = .234
  theta = 0.345
  sigma2 = .2

  # Second Derivative w.r.t phi
  d2phi = (1/((-1 + phi)^5*(1 + phi)^3*tau^2))*
    (2*sigma2*(-12*(1 + phi)^2*
                 ((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau/2) +  phi^tau) - (1/2)*(1 + theta)^2*(-1 + phi^2)*tau) + 
                 (-1 + phi)^2*(-2*(-1 + theta)^2*(3 - 4*phi^(tau/2) + phi^tau) + phi^(tau/2 - 2)*(-1 + phi^2)*(-2 + phi^(tau/2))*
                                 (-phi - theta^2*phi + theta*(1 + 4*phi + phi^2))*
                                 tau + phi^(tau/2 - 2)*(1 + phi)^2*(theta + phi + theta^2*phi + theta*phi^2)*(-1 + phi^(tau/2))*tau^2) + 6*(-1 + phi)*(1 + phi)*
                 ((theta + phi)*(1 + theta*phi)*(3 - 4*phi^(tau/2) + phi^tau) + (1/2)*(1 + theta)^2*
                    (-1 + phi^2)*tau + (1 + phi)*((-theta)*(theta + phi)*
                                                    (3 - 4*phi^(tau/2) + phi^tau) - (1 + theta*phi)*
                                                    (3 - 4*phi^(tau/2) + phi^tau) - (1 + theta)^2*phi*tau - phi^(-1 + tau/2)*(theta + phi)*
                                                    (1 + theta*phi)*(-2 + phi^(tau/2))*tau))))

  # Second derivative w.r.t theta
  d2theta = (2*sigma2*(2*phi*(3 - 4*phi^(tau/2) + phi^tau) + (-1 + phi^2)*tau))/((-1 + phi)^3*(1 + phi)*tau^2)

  # Second derivatie w.r.t sigma2
  d2sigma2 = 0

  # Partial derivative w.r.t theta and phi
  dthetadphi = (-(1/((-1 + phi)^4*(1 + phi)^2*tau^2)))*2*sigma2*(2*(3 - 4*phi^(tau/2) + phi^tau)*
                                                                   (1 + phi*(3 + phi + phi^2) + theta*(1 + phi*(2 + 3*phi))) + (2*(1 + theta)*(-1 + phi)*(1 + phi)^2 + 2*phi^(tau/2 - 1)*(-1 + phi^2)*
                                                                                                                                  (1 + 2*theta*phi + phi^2) - phi^(tau - 1)*(-1 + phi^2)*(1 + 2*theta*phi + phi^2))*tau)
  # Partial derivative w.r.t. sigma2 and phi
  dsigma2dphi = (1/((-1 + phi)^4*(1 + phi)^2*tau^2))*
    (2*((-1 + phi)*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau/2) + phi^tau) - (1/2)*(1 + theta)^2*(-1 + phi^2)*tau) +
          3*(1 + phi)*((-(theta + phi))*(1 + theta*phi)*(3 - 4*phi^(tau/2) + phi^tau) - (1/2)*(1 + theta)^2*(-1 + phi^2)*tau) - 
          (-1 + phi)*(1 + phi)*((-theta)*(theta + phi)*(3 - 4*phi^(tau/2) + phi^tau) -
                                  (1 + theta*phi)*(3 - 4*phi^(tau/2) + phi^tau) - 
                                  (1 + theta)^2*phi*tau -
                                  phi^(-1 + tau/2)*(theta + phi)*(1 + theta*phi)*(-2 + phi^(tau/2))*tau)))
  
  # Partial derivative w.r.t sigma2 and theta
  dsigma2dtheta = (2*((theta + 1)*(phi^2 - 1)*tau + (2*theta*phi + phi^2 + 1)*
                         (phi^tau - 4*phi^(tau/2) + 3)))/((phi - 1)^3*(phi + 1)*tau^2)
   
  # Test C++ is equal to deriv 2nd AR1
  expect_equal(deriv_2nd_arma11(phi, 0, sigma2, tau)[,c(1,5,3)], deriv_2nd_ar1(phi, sigma2, tau), check.attributes = F)
  
  # Test C++ is equal to deriv 2nd MA1
  expect_equal(deriv_2nd_arma11(0, theta, sigma2, tau)[,c(2,6,3)], deriv_2nd_ma1(theta, sigma2, tau), check.attributes = F)
  
  # Test C++ is equal to R implementation
  testmat = cbind(d2phi, d2theta, d2sigma2, dthetadphi, dsigma2dphi, dsigma2dtheta)

  expect_equal(deriv_2nd_arma11(phi, theta, sigma2, tau), testmat, check.attributes = F)
})