context("Convert precision matrix to partial correlation")

test_that("Conversion from precision matrix to partial correlation matrix works", {
  
  p <- 100
  xx <- diag(p)*2
  xx[1,p] <- 2
  xx[p,1] <- -2
  
  a <- Precision2PartialCor(xx)
  expect_equal(sum(a), p)
  
})