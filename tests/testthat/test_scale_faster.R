context("scale_faster")

test_that("scale and scale_faster return the same result",{
  set.seed(1)
  X  <- matrix(rnorm(20),4,5)
  Y1 <- scale(X)
  Y2 <- scale_faster(X)
  expect_equivalent(Y1,Y2,tolerance = 1e-15,scale = 1)
})
