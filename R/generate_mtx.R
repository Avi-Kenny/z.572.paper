#' Generate CAR/DAGAR model matrix components
#'
#' @param model One of: "iCAR", "CAR", "Scaled iCAR", "DAGAR", "DAGAR_OF"
#' @param adj_mtx An adjacency matrix
#' @return A list containing the following: \cr
#'     * `D`: Distance matrix for CAR-type matrices \cr
#'     * `n_i`: Vector of number of neighbors, indexed by region i \cr
#'     * `Q`: Precision matrix (with tau_W param set to 1) \cr
#'     # `neighbors`: A list of neighbors, indexed by region i
#' @export
generate_mtx <- function(model, adj_mtx, rho=1) {

  if (model %in% c("iCAR", "CAR", "Scaled iCAR")) {

    # Create neighbor sets
    neighbors <- list()
    dim <- length(adj_mtx[1,])
    n_i <- c()
    for (i in 1:dim) {
      nset <- c()
      for (j in 1:dim) {
        if (adj_mtx[i,j]==1 & i!=j) {
          nset <- c(nset, j)
          neighbors[[as.character(i)]] <- nset
        }
      }
      n_i <- c(n_i, length(nset))
    }

    # Create D
    D <- diag(n_i)

    # Create precision
    if (model %in% c("iCAR", "Scaled iCAR")) {
      Q <- (D-0.999999*adj_mtx)
    }
    if (model=="CAR") {
      Q <- (D-(rho*adj_mtx))
    }

  }

  if (model %in% c("DAGAR", "DAGAR_OF")) {

    # !!!!! TO DO: DAGAR_OF

    # Create neighbor sets
    dim <- length(adj_mtx[1,])
    neighbors <- list()
    n_i <- c()
    for (i in 1:dim) {
      nset <- c()
      for (j in 1:dim) {
        if (j<i & adj_mtx[i,j]==1) {
          nset <- c(nset, j)
        }
      }
      if (!is.null(nset)) {
        neighbors[[as.character(i)]] <- nset
        n_i <- c(n_i, length(nset))
      } else {
        neighbors[[as.character(i)]] <- 0
        n_i <- c(n_i, 0)
      }
    }

    tau <- (1+(n_i-1)*(rho^2)) / (1-(rho^2))
    FF <- diag(tau)

    B <- matrix(0, nrow=dim, ncol=dim)
    for (i in 1:dim) {
      for (j in 1:dim) {
        B[i,j] <- ifelse(
          j<i && adj_mtx[i,j]==1,
          rho / (1+((n_i[i]-1)*(rho^2))),
          0
        )
      }
    }

    L <- diag(rep(1,dim)) - B
    Q <- t(L) %*% FF %*% L
    D = NULL

  }

  return (list(
    D = D,
    n_i = n_i,
    Q = Q,
    neighbors = neighbors
  ))

}
