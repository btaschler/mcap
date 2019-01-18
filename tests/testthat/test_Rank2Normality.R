context("Rank to normality transform")

test_that("Transformation of matrix to normality based on rank works", {
  
  n <- 200
  p <- 10
  xx <- matrix(rweibull(n*p, shape = 1), n, p)
  xx_transf <- Rank2Normality(xx)
  
  for(j in seq(p)){
    a <- shapiro.test(xx_transf[,j])
    expect_equal((a$p.value>0.999), TRUE)
  }
  
})