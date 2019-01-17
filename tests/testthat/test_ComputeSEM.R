context("Compute SEM")

test_that("Computation of SEM works", {
  
  n <- 10
  x <- rnorm(n)
  a <- ComputeSEM(x)
  expect_equal(a, sd(x)/sqrt(n))
  
  n <- 20
  num_na <- 5
  x <- rnorm(n)
  x[sample(n, num_na)] <- NA
  a <- ComputeSEM(x)
  expect_equal(a, sd(x, na.rm = TRUE)/sqrt(n-num_na))
  
})