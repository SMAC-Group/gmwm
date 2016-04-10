context("Process to WV Haar - Unit Tests")

# Create a fake "zero" length vector
a = 1
a = a[0]

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
  sigma2 = .2

  expect_equal(arma11_to_wv(phi,theta,sigma2,tau), arma_to_wv(phi,theta,sigma2,tau) )
  
  # Test using ma1_to_* or ar1_to_* with arma11_to_wv()
  expect_equal(arma11_to_wv(0, theta, sigma2,tau), ma1_to_wv(theta,sigma2,tau))
  expect_equal(arma11_to_wv(phi, 0, sigma2,tau), ar1_to_wv(phi,sigma2,tau))

})