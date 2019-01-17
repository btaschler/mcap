context("Centre data per group")

test_that("Group-specific mean-centering of data matrix works", {
  
  n <- 200
  p <- 100
  num_groups <- 2
  xx <- matrix(rnorm(n*p), n, p)
  labels <- ceiling(runif(n) * num_groups)
  
  a <- CentrePerGroup(xx, true_labels = labels)
  for(k in seq(num_groups)){
    expect_equal(c(a[labels==k, ]), 
                 c(scale(xx[labels==k, ], center = TRUE, scale = FALSE)))
  }
  
})