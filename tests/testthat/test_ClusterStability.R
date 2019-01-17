context("Cluster stability")

test_that("Computation of cluster stability works", {
  
  # Easy clustering task with two very distinct groups. 
  # Info: The test is based on the expectation that the adjusted Rand index is 1 for 
  # every pair-wise clustering on two subsets of the data.
  
  n <- 100
  p <- 2
  num_groups <- 2
  labels <- floor(runif(n) * num_groups)
  xx <- matrix(rnorm(n*p) + 10*labels, n, p)
  a <- ClusterStability(xx, k = num_groups, B = 5, frac_subsample = 0.75)
  expect_equal(sum(a), length(a))
  
})