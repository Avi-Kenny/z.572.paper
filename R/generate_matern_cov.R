#' Generate a matern covariance matrix from a graph
#'
#' @param type Type of dataset (either "path", "grid", or "US")
#' @param phi Covariance matrix parameter
#' @return A covariance matrix
#' @export
generate_matern_cov <- function(type, phi) {

  n <- ifelse(type=="US", 48, 100)
  matern <- matrix(NA, nrow=n, ncol=n)

  if (type=="path") {

    for (i in 1:n) {
      for (j in 1:n) {
        matern[i,j] <- exp(-1*phi*abs(i-j))
      }
    }

  }

  if (type=="grid") {

    coords <- list( x=rep(1:10, 10), y=rep(1:10, each=10) )

    for (i in 1:n) {
      for (j in 1:n) {
        dist <- sqrt((coords$x[i]-coords$x[j])^2+(coords$y[i]-coords$y[j])^2)
        matern[i,j] <- exp(-1*phi*dist)
      }
    }

  }

  if (type=="US") {

    adj_mtx_us <- adj_from_gis_us()
    m <- adj_mtx_us * mtx_dist_US
    num_neighbors <- sum(adj_mtx_us)

    avg_nbr_dist <- sum(m) / num_neighbors

    scaled_dist_mtx <- mtx_dist_US / avg_nbr_dist

    for (i in 1:n) {
      for (j in 1:n) {
        dist <- scaled_dist_mtx[i,j]
        matern[i,j] <- exp(-1*phi*dist)
      }
    }

  }

  return (matern)

}
