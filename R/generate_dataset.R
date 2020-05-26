#' Generate one dataset !!!!! TO DO
#'
#' @param type Type of dataset (either "path", "grid", or "US")
#' @param tau_w Spatial parameter
#' @param cov_type Type of covariance matrix ("Matern", "DAGAR", or "CAR")
#' @param phi Matern covariance matrix parameter
#' @param rho DAGAR/CAR covariance matrix parameter
#' @return A list containing w, y, x1, x2
#' @export
generate_dataset <- function(type, tau_w, cov_type="Matern", phi=NA, rho=NA) {

  k <- ifelse(type=="US", 48, 100)

  # Hard-code parameters
  beta1 <- 1
  beta2 <- 5
  tau_e <- 2.5

  # Sample covariates and residual error
  x1 <- rnorm(n=k)
  x2 <- rnorm(n=k)
  e <- rnorm(n=k, mean=0, sd=1/sqrt(tau_e))

  # Generate covariance matrix
  if (cov_type=="Matern") {

    M_inv <- generate_matern_cov(type=type, phi=phi)
    cov_mtx <- (1/tau_w)*M_inv

  } else {

    if (type=="path") { adj_mtx <- generate_graph_path() %>% adj_from_graph() }
    if (type=="grid") { adj_mtx <- generate_graph_grid() %>% adj_from_graph() }
    if (type=="US") { adj_mtx <- adj_from_gis_us() }

  }

  # Generate DAGAR matrix
  if (cov_type=="DAGAR") {

    mtx <- generate_mtx(model="DAGAR", adj_mtx=adj_mtx, rho=rho)
    cov_mtx <- (1/tau_w)*solve(mtx$Q)

    }

  # Generate CAR matrix
  if (cov_type=="CAR") {

    mtx <- generate_mtx(model="CAR", adj_mtx=adj_mtx, rho=rho)
    cov_mtx <- (1/tau_w)*solve(mtx$Q)

  }
  w <- as.numeric(rmvnorm(n=1, sigma=cov_mtx))

  # Generate data vector
  y <- beta1*x1 + beta2*x2 + w + e

  return (list(
    y = y,
    w = w,
    x1 = x1,
    x2 = x2
  ))

}
