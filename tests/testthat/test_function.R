library(testthat)        # load testthat package
library(tnlTEST)            # load our package

# Test whether the output is a list
test_that("functions returns a list", {
  x=c(1.0021,-0.0128,2.739,1.605,1.084,3.072,2.417)
  y=c(-1.159,-0.038,-0.374,0.900,-1.1063,2.739,3.104)
  expect_type(tnl(7,1), "list")
  expect_type(tnl.sim(5,1), "list")
  expect_type(ptnl(2,5,1), "list")
  expect_type(ptnl(2,5,1,exact="TRUE"), "list")
  expect_type(ptnl(2,5,1,exact="FALSE"), "list")
  expect_type(ptnl(2,11,1), "list")
  expect_type(dtnl(2,5,1), "list")
  expect_type(dtnl(2,5,1,exact="TRUE"), "list")
  expect_type(dtnl(2,5,1,exact="FALSE"), "list")
  expect_type(dtnl(2,11,1), "list")
  expect_type(qtnl(.2,5,1), "list")
  expect_type(qtnl(.2,5,1,exact="TRUE"), "list")
  expect_type(qtnl(.2,5,1,exact="FALSE"), "list")
  expect_type(qtnl(.2,11,1), "list")
  expect_type(tnl.test(x,y,l=2), "list")
  expect_type(tnl.test(x,y,l=2,exact="TRUE"), "list")
  expect_type(tnl.test(x,y,l=1,exact="FALSE"), "list")
  expect_type(tnl.test(c(x,0.975,1.144,0.572,-0.532),c(y, 1.007,-2.023,1.468,1.396),l=2), "list")
})
## Test whether the output return the right number
test_that("functions returns the right output", {
  expect_equal(sum(tnl(7,1)$pmf), 1)
  expect_equal(sum(tnl.sim(4,1)$pmf), 1)
  expect_equal(tnl(7,1)$cdf[7], 1)
  expect_equal(tnl.sim(4,1)$cdf[4], 1)
  expect_equal(ptnl(5,5,2)$cdf, 1)
})
## Test whether the output contains the right number
test_that("functions returns a list with the specified length", {
  expect_length(tnl(7,1)$pmf, 7)
  expect_length(tnl(7,1)$cdf, 7)
  expect_length(tnl.sim(5,1)$pmf, 5)
  expect_length(tnl.sim(5,1)$cdf, 5)
  expect_length(ptnl(2,7,1)$cdf, 1)
  expect_length(dtnl(2,7,1)$pmf, 1)
  expect_length(ptnl(2,11,1)$cdf, 1)
  expect_length(dtnl(2,11,1)$pmf, 1)
  expect_length(qtnl(.2,7,1)$quantile, 1)
})
## Test whether the output is a vector with the expected size
test_that("functions returns a  vector with the expected size", {
  expect_vector(rtnl(10,7,2), ptype = double(), size = 10)
  expect_vector(dtnl.lehmann(2,6,0.8), ptype = double(), size = 6)
  expect_vector(ptnl.lehmann(2,5,1.2), ptype = double(), size = 5)
})
## Test whether the output should not exceed one.
test_that("functions returns number should not exceed one", {
  x=c(1.0021,-0.0128,2.739,1.605,1.084,3.072,2.417)
  y=c(-1.159,-0.038,-0.374,0.900,-1.1063,2.739,3.104)
  expect_lte(tnl(7,1)$pmf[2], 1)
  expect_lte(tnl(7,1)$cdf[5], 1)
  expect_lte(tnl.sim(5,1)$pmf[4], 1)
  expect_lte(tnl.sim(5,1)$cdf[3], 1)
  expect_lte(ptnl(2,7,1)$cdf, 1)
  expect_lte(dtnl(2,7,1)$pmf, 1)
  expect_lte(tnl.test(x,y,l=3)$p.value, 1)
})
## Test whether the output should exceed zero.
test_that("functions returns number should exceed zero", {
  x=c(1.0021,-0.0128,2.739,1.605,1.084,3.072,2.417)
  y=c(-1.159,-0.038,-0.374,0.900,-1.1063,2.739,3.104)
  expect_lte(0, tnl(7,1)$pmf[5])
  expect_lte(0, tnl(7,1)$cdf[2])
  expect_lte(0, tnl.sim(5,1)$pmf[4])
  expect_lte(0,tnl.sim(5,1)$cdf[3])
  expect_lte(0, ptnl(2,7,1)$cdf)
  expect_lte(0,dtnl(2,7,1)$pmf)
  expect_lte(0,tnl.test(x,y,l=3)$p.value)
})

## Test whether the code throw an error.
test_that("functions returns errors", {
  x=c(1.0021,-0.0128,2.739,1.605,1.084,3.072,2.417)
  y=c(-1.159,-0.038,-0.374,0.900,-1.1063,2.739,3.104)
  expect_error(tnl(4,2),"n must be > 2l")
  expect_error(tnl.sim(4,2),"n must be > 2l")
  expect_error(ptnl(1,7,2),"q must be >= l" )
  expect_error(ptnl(3,4,2),"n must be > 2l")
  expect_error(dtnl(1,7,2),"k must be >= l" )
  expect_error(dtnl(3,4,2),"n must be > 2l")
  expect_error(qtnl(.3,4,2),"n must be > 2l")
  expect_error(qtnl(1.3,8,1),"p must be between 0 and 1")
  expect_error(dtnl.lehmann(3,5,1.2),"n must be > 2l")
  expect_error(ptnl.lehmann(3,5,1.2),"n must be > 2l")
  expect_error(tnl.test(x,c(y,2.1),l=2),"The length of x and y must be equal")
  expect_warning(tnl.test(c(x,NA),c(NA,y),l=2),"Since the data should not contain missing values, we exclude the missing values from the data")
})
