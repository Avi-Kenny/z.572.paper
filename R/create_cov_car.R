#' Create a covariange matrix (CAR model) from an adjacency matrix
#'
#' @param adj_mtx Adjacency matrix
#' @param rho Spatial correlation parameter
#' @return A covariance matrix
#' @export
generate_cov_car <- function(adj_mtx, rho) {

  dim <- length(adj_mtx[1,])

  # Create neighbor sets
  neighbors <- list()
  n_lti <- c()
  for (i in 1:dim) {
    nset <- c()
    for (j in 1:dim) {
      if (adj_mtx[i,j]==1) {
        nset <- c(nset, j)
      }
    }
    neighbors[[as.character(i)]] <- nset
    n_lti <- c(n_lti, length(nset))
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
