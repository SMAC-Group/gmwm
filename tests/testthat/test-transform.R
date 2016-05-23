context("Parameter Transformations - Unit Tests")

# Create a fake "zero" length vector
a = numeric()

# Overall taus
tau = 2^(1:4)

test_that("Transformation - Logit",{
  
  # Check Version Conversion
  v = as.matrix( c(.1, .9,.3) )
  
  a = logit(v)
  
  b = logit_inv(a)
  
  expect_equal(b,v,check.attributes = F)
  
  # Check Inf handling
  
  v = as.matrix( c(-.1, 1.1, 1.2) )
  
  a = logit(v)
  
  expect_true(all(is.nan(a)))
  
  # No Recovery
  b = logit_inv(a)
  
  expect_true(all(is.nan(b)))
})


test_that("Transformation - Pseudo Logit",{
  
  # Check Version Conversion
  v = as.matrix( c(.1, .9,.3) )
  
  a = pseudo_logit(v)
  
  b = pseudo_logit_inv(a)
  
  expect_equal(b,v,check.attributes = F)
  
  # Check Inf handling
  v = as.matrix( c(-2, -1.1, 1.1) )
  
  a = pseudo_logit(v)
  
  expect_true(all(is.nan(a)))
  
  # No Recovery
  b = pseudo_logit_inv(a)
  
  expect_true(all(is.nan(b)))
})

test_that("Transformation - logit2",{
  
  # Check Version Conversion
  v = as.matrix( c(.1, .9,.3) )
  
  a = logit2(v)
  
  b = logit2_inv(a)
  
  expect_equal(b,v,check.attributes = F)
  
  # Check Inf handling
  v = as.matrix( c(-2, 0, 2) )
  
  a = logit2(v)
  
  expect_false(all(is.nan(a)))
  
  # Add Recovery
  b = logit2_inv(a)
  
  expect_equal(b,v)
  
  # Check Inf handling
  v = as.matrix( c(-3, -2.5, 2.5, 3) )
  
  a = logit2(v)
  
  expect_true(all(is.nan(a)))
  
  # Add Recovery
  b = logit2_inv(a)
  
  expect_true(all(is.nan(b)))
})



test_that("Transform Parameters",{
  
  # Settings
  phi = .62
  phi_sigma2 = .9
  wn = 2
  dr = 0.0001
  qn = 0.2
  rw = 0.75
  theta = 0.1
  theta_sigma2 = 0.5
  
  phi2 = 0.78
  theta2 = 0.32
  sigma2 = 1.25
  
  phi2arma = 0.21
  theta2arma = 0.84
  sigma2arma = 1.5
  
  # Test Overalls
  model = AR1(phi, phi_sigma2)  + MA1(theta, theta_sigma2) + 
    ARMA11(phi2, theta2, sigma2) + 
    ARMA(phi2arma,theta2arma,sigma2arma) +
    WN(wn) + DR(dr) + QN(qn) + RW(rw) 
  
  
  # Transform Values underneath the IMU Case
  cpp_transform_imu = transform_values(model$theta, model$desc, model$obj.desc, "imu")
  
  # Transform Values underneath the SSM Case
  cpp_transform_ssm = transform_values(model$theta, model$desc, model$obj.desc, "ssm")
  
  # IMU Case transformation reversion
  cpp_untransform_imu = untransform_values(cpp_transform_imu, model$desc, model$obj.desc, "imu")
  
  # SSM Case transformation reversion
  cpp_untransform_ssm = untransform_values(cpp_transform_ssm, model$desc, model$obj.desc, "ssm")
  
  
  # Check that untransform imu matches with the model statement
  expect_equal(cpp_untransform_imu, as.matrix(model$theta))
  
  # Check that untransform ssm matches with the model statement
  expect_equal(cpp_untransform_ssm, as.matrix(model$theta))
})
