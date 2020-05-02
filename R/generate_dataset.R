#' Generate one dataset !!!!! TO DO
#'
#' @param type Type of dataset (either "path", "grid", or "US")
#' @param rho Spatial parameter
#' @param tau_w Spatial parameter
#' @param phi Matern covariance matrix parameter
#' @return !!!!! A list containing the following: \cr
#'     * `params`: !!!!! A list of the parameters supplied in the function call \cr
#'     * `reg_data`: !!!!! A data frame of data points corresponding to a single
#'         region \cr
#'     * `adj_mtx`: !!!!! The adjacency matrix \cr
#' @export
generate_dataset <- function(type, rho, tau_w, phi) {

  # # !!!!! testing
  # generate_dataset(type="path", rho=0.5, tau_w=0.25, phi=0.5)

  k <- ifelse(type=="US", 48, 100)

  # Hard-code parameters
  beta1 <- 1
  beta2 <- 5
  tau_e <- 2.5

  x1 <- rnorm(n=k)
  x2 <- rnorm(n=k)
  e <- rnorm(n=k, mean=0, sd=sqrt(1/tau_e))

  w <- 999

  y <- beta1*x1 + beta2*x2 + w + e

  switch(
    type,
    "path" = {

      # !!!!! 100 vertices
      i <- 1:100
      w_i <- rep(999,100)
      n_i <- rep(999,100)

      # Construct covariance matrix (DAGAR)


    },
    "grid" = {
      # !!!!! 10 x 10 grid
    },
    "US" = {

      # !!!!!

    }
  )

  return(list(
    "params" = as.list(match.call()),
    "reg_data" = reg_data,
    "adj_mtx" = 999
  ))

}
