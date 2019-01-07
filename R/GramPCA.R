
GramPCA <- function(xx, npc) {
  #' Perform PCA using Gram matrix
  #'
  #' This function performs PCA using the eigenvectors of the Gram matrix G 
  #' (G=X*X') when n < p, or standard PCA when n > p.
  #' 
  #' @author Bernd Taschler: \email{bernd.taschler@dzne.de}
  #' @author Sach Mukherjee: \email{sach.mukherjee@dzne.de}
  #' @seealso \code{\link{pcaMethods::pca}}
  #' @seealso \code{\link{eigen}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param npc The number of principal components to be returned. 
  #' 
  #' @return \item{zz}{Principal components (n x npc).}
  #' @export
  
  
  ## preliminaries
  n <- nrow(xx)
  p <- ncol(xx)
  
  ## mean centre columns
  xx <- scale(xx, center=TRUE, scale=FALSE)

  ## perform Gram-PCA or standard PCA
  if(n < p){
    gg <- xx %*% t(xx)             #compute Gram matrix
    specDecomp <- eigen(gg,npc)    #note: eigenvalues are already sorted
    lam <- specDecomp$values
    V <- specDecomp$vectors
    zz <- V[,1:npc]
    
  }else{ #case n > p
    pcaS4 <- pcaMethods::pca(xx, nPcs=npc, method='svd', scale='none')
    zz <- pcaS4@scores 
  }
  
  return(list('zz' = zz))
}