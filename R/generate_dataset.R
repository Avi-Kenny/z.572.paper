#' Generate one dataset !!!!! TO DO
#'
#' @param type Type of dataset (either "path", "grid", or "US")
#' @param tau_w Spatial parameter
#' @param phi Matern covariance matrix parameter
#' @return !!!!! A list containing the following: \cr
#'     * `params`: !!!!! A list of the parameters supplied in the function call \cr
#'     * `reg_data`: !!!!! A data frame of data points corresponding to a single
#'         region \cr
#'     * `adj_mtx`: !!!!! The adjacency matrix \cr
#' @export
generate_dataset <- function(type, tau_w, phi) {

  k <- ifelse(type=="US", 48, 100)

  # Hard-code parameters
  beta1 <- 1
  beta2 <- 5
  tau_e <- 2.5

  # Sample covariates and residual error
  # !!!!! May want to pre-calculate M_inv
  x1 <- rnorm(n=k)
  x2 <- rnorm(n=k)
  e <- rnorm(n=k, mean=0, sd=1/sqrt(tau_e))
  M_inv <- generate_matern_cov(type=type, phi=phi)
  w <- as.numeric(rmvnorm(n=1, sigma=(1/tau_w)*M_inv))

  # Generate data vector
  y <- beta1*x1 + beta2*x2 + w + e

  return (list(
    y = y,
    w = w,
    x1 = x1,
    x2 = x2
  ))

}
