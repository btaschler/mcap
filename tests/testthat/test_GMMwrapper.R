context("Wrapper for GMM")

test_that("Running nethet::mixglasso() works", {
  
  # Easy clustering scenario: expect aRI=1
  n <- 100
  p <- 10
  num_groups <- 2
  labels <- floor(runif(n) * num_groups)
  xx <- matrix(rnorm(n*p) + 10*labels, n, p)
  a <- GMMwrapper(xx, k = num_groups, true_labels = labels)
  expect_equal(a$aRI, 1)
  
})