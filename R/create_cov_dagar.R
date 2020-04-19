#' Create a covariange matrix (DAGAR model) from an adjacency matrix
#'
#' @param adj_mtx Adjacency matrix
#' @param rho Spatial correlation parameter
#' @return A list containing the following: \cr
#'     * `mtx_cov`: Covariance matrix corresponding to the given adjacency matrix
#'     * `mtx_B`: B matrix
#'     * `n_lti`: Vector containing number of neighbors for region i
#'     * `neighbors`: List containing neighbor sets for each region
#'     * `tau`: Vector of tau_1, ..., tau_k
#'     * `dim`: Dimension (number of regions)
#' @export
generate_cov_dagar <- function(adj_mtx, rho) {

  dim <- length(adj_mtx[1,])

  # Create neighbor sets
  neighbors <- list()
  n_lti <- c()
  for (i in 1:dim) {
    nset <- c()
    for (j in 1:dim) {
      if (j<i & adj_mtx[i,j]==1) {
        nset <- c(nset, j)
      }
    }
    if (!is.null(nset)) {
      neighbors[[as.character(i)]] <- nset
      n_lti <- c(n_lti, length(nset))
    } else {
      neighbors[[as.character(i)]] <- 0
      n_lti <- c(n_lti, 0)
    }
  }

  # Create F matrix
  # !!!!! Make sure tau_1 is supposed to equal 1
  tau <- (1 + (n_lti-1)*(rho^2)) / (1-rho^2)
  mtx_F <- diag(tau)

  # Create B matrix
  mtx_B <- matrix(0, nrow=dim, ncol=dim)
  for (i in 2:dim) {
    for (j in 1:dim) {
      if (j %in% neighbors[[as.character(i)]]) {
        mtx_B[i,j] <- rho / (1 + (n_lti[i]-1)*(rho^2))
      }
    }
  }

  # Create covariance matrix
  mtx_I <- diag(rep(1,dim))
  mtx_L <- mtx_I - mtx_B
  mtx_cov <- solve( t(mtx_L) %*% mtx_F %*% mtx_L )

  return (list(
    "mtx_cov" = mtx_cov,
    "mtx_B" = mtx_B,
    "n_lti" = n_lti,
    "neighbors" = neighbors,
    "tau" = tau,
    "dim" = dim
  ))

}
