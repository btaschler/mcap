#' mcap: Model-based clustering in very high dimensions via adaptive projections.
#' 
#' Model-based clustering in very high dimensions (especially p >> n) via adaptive projections. 
#' Clustering is based on Gaussian mixture models in a lower dimensional (projected) space.
#' Projection dimension is set adaptively based on a cluster stability criterion.
#' Available projection variants (so far) include PCA and Random Projections.
#'
#' @section Main mcap functions:
#'   \code{"MCAPfit"}
#'   \code{"GMMwrapper"}
#'   \code{"OptDimClusterStability"}
#'   \code{"ClusterStability"}
#'   \code{"GramPCA"}
#'
#' @name mcap
#' @docType package
NULL
