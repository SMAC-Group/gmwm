context("Graphical Functions for WV, AV and HV - Unit Tests")

## 'avar'

# 'avar' data structure:
# Update autoplot.avar if this test fails
test_that("avar data structure",{
  # create a 'avar' object
  set.seed(999)
  x=rnorm(100)
  av = avar(x)
  av.name = names(av)
  exp.name = c("clusters", "allan", "errors", "adev", "lci", "uci", "type"  )
  
  expect_equal(av.name, exp.name)
  
})


## 'hadam'

# 'hadam' data structure:
# Update autoplot.hadam if this test fails
test_that("hadam data structure",{
  # create a 'hadam' object
  set.seed(999)
  x=rnorm(100)
  hv = hadam(x)
  hv.name = names(hv)
  exp.name = c("clusters", "hadamard", "errors", "hdev", "lci", "uci", "type")
  
  expect_equal(hv.name, exp.name)
  
})
