
RandProject <- function(xx, q, method = 'gaussian'){
  #' Random projection of input matrix X
  #'
  #' Perform a random projection (Gauss, Achlioptas or Li) of the input matrix X.
  #' Based on the RandPro package. 
  #' 
  #' @seealso \code{\link{RandPro::form_matrix}}
  #' 
  #' @param xx The data matrix (n x p).
  #' @param q The target dimension. 
  #' @param method Projection method using either a Gaussian (\code{"gaussian"}), 
  #'               a sparse (\code{"achlioptas"}) or very sparse (\code{"li"}) projection 
  #'               matrix. Default: \code{"gaussian"}.
  #' 
  #' @return Projected matrix (n x q).
  #' @export
  
  p <- ncol(xx)
  rr <- RandPro::form_matrix(rows = p, cols = q, JLT = FALSE, 
                             eps = 0.1, projection =  method)

  return(xx %*% rr)
}