#' Converts Fahrenheit to Kelvin
#'
#' This function performs PCA using the eigenvectors of the Gram matrix G 
#' (G=X*X') when n < p, or standard PCA when n > p.
#' @param xx The data matrix (n x p).
#' @param npc The number of principal components to be returned. 
#' @return Principal components (n x npc).
#' @export
#' @examples

GramPCA <- function(xx, npc) {
  ## preliminaries
  n <- nrow(xx)
  p <- ncol(xx)
  
  ## mean centre columns
  xx <- scale(xx, center=TRUE, scale=FALSE)

  ## compute Gram matrix
  gg <- xx %*% t(xx)
  
  ## perform Gram-PCA or standard PCA
  if(n < p){
    specDecomp <- eigen(gg,npc)
    lam <- specDecomp$values
    V <- specDecomp$vectors
    #lamSorted <- sort.int(lam, decreasing = T, index.return = T) #eigenvals are already sorted
    zz <- V[,1:npc]
    
  }else{ #case n > p
    pcaS4 <- pcaMethods::pca(xx, nPcs=npc, method='svd', scale='none')
    zz <- pcaS4@scores 
  }
  
  return(zz)
}