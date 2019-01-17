context("Random projections")

test_that("Random projection of matrix works", {
  
  n <- 200
  p <- 100
  q <- 10
  xx <- matrix(rnorm(n*p), n, p)
  a <- RandProject(xx, q, method = 'gaussian')
  expect_equal(dim(a), c(n,q))
  
  a <- RandProject(xx, 2*q, method = 'achlioptas')
  expect_equal(dim(a), c(n,2*q))
  
})