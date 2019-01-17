context("Run MCAP")

test_that("Full MCAP model fit works", {
  
  # Easy clustering scenario: expect aRI=1
  n <- 100
  p <- 10
  num_groups <- 2
  labels <- floor(runif(n) * num_groups)
  xx <- matrix(rnorm(n*p) + 10*labels, n, p)
  
  a <- MCAPfit(xx, k = num_groups, projection = 'PCA', true_labels = labels)
  expect_equal(a$fit_q_opt$q_opt, 2)
  expect_equal(a$fit_gmm$aRI, 1)
  
  a <- MCAPfit(xx, k = num_groups, projection = 'achlioptas', true_labels = labels)
  expect_equal(a$fit_q_opt$q_opt, 2)
  expect_equal(a$fit_gmm$aRI, 1)
  
})