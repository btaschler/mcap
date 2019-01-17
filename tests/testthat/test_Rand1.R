context("Random sign vector")

test_that("Random sign vector has only +1/-1 elements", {
  
  n <- 1000
  s <- Rand1(n)
  mean(s)
  
  expect_equal(round(mean(s)), 0)
  expect_equal(mean(sort(unique(s)) == c(-1,1)), 1)
  
})