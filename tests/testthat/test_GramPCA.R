context("Gram PCA")

test_that("PCA works if n < p", {
  
  xx <- scale(matrix(rnorm(200),10,20), center=TRUE, scale=FALSE)
  a <- GramPCA(xx, npc = 1)$zz
  b <- pcaMethods::pca(xx, nPcs=1, method='svd', scale='none')@scores
  expect_equal(a, c(b * a[1]/b[1]))
  
  xx <- scale(matrix(runif(2000)*100,10,200), center=TRUE, scale=FALSE)
  a <- GramPCA(xx, npc = 1)$zz
  b <- pcaMethods::pca(xx, nPcs=1, method='svd', scale='none')@scores
  expect_equal(a, c(b * a[1]/b[1]))
})