#' Create a covariange matrix (CAR model) from an adjacency matrix
#'
#' @param adj_mtx Adjacency matrix
#' @param rho Spatial correlation parameter
#' @return A list containing the following: \cr
#'     * `mtx_cov`: Covariance matrix corresponding to the given adjacency
#'         matrix
#'     * `mtx_D`: D matrix
#'     * `n_i`: Vector containing number of neighbors for region i
#'     * `neighbors`: List containing neighbor sets for each region
#'     * `dim`: Dimension (number of regions)
#' @export
generate_cov_car <- function(adj_mtx, rho) {

  dim <- length(adj_mtx[1,])

  # Create neighbor sets
  neighbors <- list()
  n_i <- c()
  for (i in 1:dim) {
    nset <- c()
    for (j in 1:dim) {
      if (adj_mtx[i,j]==1 & i!=j) {
        nset <- c(nset, j)
      }
    }
    neighbors[[as.character(i)]] <- nset
    n_i <- c(n_i, length(nset))
  }

  # Create covariance matrix
  mtx_D <- diag(n_i)
  mtx_cov <- solve( mtx_D - rho*adj_mtx )

  return (list(
    "mtx_cov" = mtx_cov,
    "mtx_D" = mtx_D,
    "n_i" = n_i,
    "neighbors" = neighbors,
    "dim" = dim
  ))

}
