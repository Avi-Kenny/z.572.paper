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
generate_mtx <- function(model, adj_mtx, rho=0.99, p_order=NA) {

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

  if (model == "DAGAR") {

    # If necessary, permute adjacency matrix
    if (!is.na(p_order)) {
      mtx_p <- as(as.integer(p_order), "pMatrix")
      adj_mtx <- mtx_p %*% adj_mtx %*% t(mtx_p)
    }

    # Get number of neighbors
    dim <- length(adj_mtx[1,])
    neighbors <- list()
    mtx_nbrs <- adj_mtx * lower.tri(adj_mtx)
    n_i <- as.numeric(mtx_nbrs %*% rep(1,dim))

    # Get neighbors
    mtx_cols <- matrix(rep(1:dim, each=dim), nrow=dim)
    mtx_nbrs2 <- mtx_nbrs * mtx_cols
    for (i in 1:dim) {
      neighbors[[as.character(i)]] <- mtx_nbrs2[i,][mtx_nbrs2[i,]!=0]
    }

    # Create matrices
    tau <- (1+(n_i-1)*(rho^2)) / (1-(rho^2))
    FF <- diag(tau)
    b_i <- rho / (1+((n_i-1)*(rho^2)))
    B <- mtx_nbrs * matrix(rep(b_i,dim), nrow=dim)
    L <- diag(rep(1,dim)) - B
    Q <- t(L) %*% FF %*% L
    D <- NULL

    # # Permute back vectors
    # if (!is.na(p_order)) {
    #   neighbors_new <- list()
    #   for (i in 1:length(neighbors)) {
    #     neighbors_new[[i]] <- neighbors[[p_order[i]]]
    #   }
    #   neighbors <- neighbors_new
    #   n_i <- t(mtx_p) %*% n_i
    # }

  }

  if (model == "DAGAR_OF") {

    dim <- length(adj_mtx[1,])
    Q_list <- list()
    counter <- 0

    # DAGAR_OF: this is a stochastic approximation
    for (p in list(1:dim, dim:1, sample(1:dim), sample(1:dim))) {

      # Increment counter
      counter <- counter+1

      # Permute adjacency matrix
      mtx_p <- as(as.integer(p), "pMatrix")
      mtx_p_inv <- Matrix::t(mtx_p)
      adj_mtx_p <- mtx_p %*% adj_mtx %*% mtx_p_inv

      # Get number of neighbors
      # neighbors <- list()
      mtx_nbrs <- adj_mtx_p * lower.tri(adj_mtx_p)
      n_i <- as.numeric(mtx_nbrs %*% rep(1,dim))

      # Get neighbors
      # mtx_cols <- matrix(rep(1:dim, each=dim), nrow=dim)
      # mtx_nbrs2 <- mtx_nbrs * mtx_cols
      # for (i in 1:dim) {
      #   neighbors[[as.character(i)]] <- mtx_nbrs2[i,][mtx_nbrs2[i,]!=0]
      # }

      # Create matrices
      tau <- (1+(n_i-1)*(rho^2)) / (1-(rho^2))
      FF <- diag(tau)
      b_i <- rho / (1+((n_i-1)*(rho^2)))
      B <- mtx_nbrs * matrix(rep(b_i,dim), nrow=dim)
      L <- diag(rep(1,dim)) - B
      Q_list[[counter]] <- mtx_p_inv %*% t(L) %*% FF %*% L %*% mtx_p

    }

    n_i <- NULL
    neighbors <- NULL
    D <- NULL
    Q <- (Q_list[[1]]+Q_list[[2]]+Q_list[[3]]+Q_list[[4]])/4

  }

  return (list(
    D = D,
    n_i = n_i,
    Q = Q,
    neighbors = neighbors
  ))

}
