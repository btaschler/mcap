
MCAPfit <- function(xx, k, projection = 'PCA',
                    true_labels = NULL, centering_per_group = FALSE, 
                    parallel = FALSE, verbose = FALSE, ...){
  #' Model based clustering via adaptive (linear) projections
  #'
  #' \code{"MCAPfit"} performs model based clustering using full variance 
  #' Gaussian mixtures in a lower dimensional projected space obtained via 
  #' adaptive (linear) projections.
  #' Projection variants include PCA-based and random projection. 
  #' 
  #' @author Bernd Taschler: \email{bernd.taschler@dzne.de}
  #' @author Sach Mukherjee: \email{sach.mukherjee@dzne.de}
  #' @references Taschler, B., Dondelinger, F. and Mukherjee, S. (2019) 
  #'             Model based clustering via adaptive projections \url{https://arxiv.org/pdf/??.pdf}
  #' @seealso \code{\link{GMMwrapper}}, 
  #'          \code{\link{OptDimClusterStability}},
  #'          \code{\link{ClusterStability}} 
  #' 
  #' @param xx The data matrix (n x p).
  #' @param k The number of clusters.
  #' @param projection Projection method (\code{"PCA"}, \code{"gaussian"}, 
  #'                   \code{"achlioptas"} or \code{"li"}). Default: \code{"PCA"}.
  #' @param true_labels Vector of true cluster assignments (when provided, it is 
  #'                    used to compute the Rand index). 
  #' @param centering_per_group Logical, when true: mean centre input matrix (if true 
  #'                            labels are provided: centre data per group)
  #' @param parallel Logical, when true: perform line search over projection 
  #'                 dimension in parallel. 
  #' @param verbose Logical, when true: print some progress information. 
  #' @param ... Additional options for \code{\link{OptDimClusterStability}} and
  #'            \code{\link{GMMwrapper}}. 
  #' 
  #' @return \item{fit_gmm}{ Model fit (GMM output of \code{\link[nethet]{mixglasso}}), 
  #'                         including BIC, MMDL and adj. Rand index (when 
  #'                         \code{true_labels} is provided).}
  #'         \item{fit_q_opt}{ Output of \code{\link{OptDimClusterStability}}.}
  #' 
  #' @details ...
  #' 
  #' @examples 
  #'   \dontrun{
  #'   MCAPfit(matrix(rnorm(200),20,10), 2, 'PCA', round(runif(20)), F, F, T)
  #'   }
  #'   
  #'   \dontrun{
  #'   MCAPfit(matrix(rnorm(5e4),100,500), 2, 'PCA', round(runif(100)), F, F, T)  #' 
  #'   }
  #' 
  #' @export
  
  ## input checks
  if(!is.null(true_labels)){ 
    stopifnot(length(true_labels) == nrow(xx)) 
    stopifnot(length(unique(true_labels)) == k)
  }
  
  ## preliminaries
  xx <- as.matrix(xx)
  n <- nrow(xx)
  p <- ncol(xx)
  
  if(centering_per_group){
    xx <- CentrePerGroup(xx, true_labels = true_labels)
  }
  
  ## determine optimal projection dimension
  fit_q_opt <- OptDimClusterStability(xx, k = k, method = projection, 
                                      true_labels = true_labels, 
                                      verbose = verbose, ...)
  
  ## parameter estimation and GMM clustering with optimised target dimension
  fit_gmm <- GMMwrapper(GramPCA(xx, npc = fit_q_opt$q_opt)$zz, 
                        k = k, true_labels = true_labels, 
                        verbose = verbose, ...)
  
  return(list('fit_gmm' = fit_gmm, 'fit_q_opt' = fit_q_opt))
}