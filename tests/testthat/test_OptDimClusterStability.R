context("Optimal projection dimension (cluster stability)")

test_that("Finding the optimal projection dimension works", {
  
  # Easy clustering scenario: expect aRI=1, q_opt=2
  n <- 100
  p <- 10
  num_groups <- 2
  labels <- floor(runif(n) * num_groups)
  xx <- matrix(rnorm(n*p) + 10*labels, n, p)
  a <- OptDimClusterStability(xx, k = num_groups, method = 'PCA', 
                              true_labels = labels)
  expect_equal(a$stab_score, 1)
  expect_equal(a$q_opt, 2)
  
})