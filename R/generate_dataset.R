#' Generate one dataset !!!!! TO DO
#'
#' @param type Type of dataset (either "path", "grid", or "US")
#' @param rho Spatial correlation parameter
#' @return A list containing the following: \cr
#'     * `params`: A list of the parameters supplied in the function call \cr
#'     * `reg_data`: A data frame of data points corresponding to a single
#'         region \cr
#'     * `adj_mtx`: The adjacency matrix \cr
#'     * `x`: X \cr
#'     * `x`: X \cr
#'     * `x`: X \cr
#'     * `x`: X \cr
#'     * `x`: X
#' @export
generate_dataset <- function(type, rho) {

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

      # Read in data
      geo <- readOGR(
        paste0(
          system.file(package = "z.572.paper"),
          "/extdata/GADM/gadm36_USA_shp"
        ),
        layer = "gadm36_USA_1",
        verbose = F
      )
      n_features <- length(geo[["GID_1"]])

      # Generate adjacency matrix
      # !!!!! Save this object for quick loading
      neighbor_list <- poly2nb(
        geo,
        queen = FALSE,
        row.names = 1:n_features
      )
      adj_mtx <- nb2mat(neighbor_list, style="B",zero.policy=TRUE)
      adj_mtx <- as.matrix(adj_mtx[1:dim(adj_mtx)[1], 1:dim(adj_mtx)[1]])


    }
  )

  return(list(
    "params" = as.list(match.call()),
    "reg_data" = reg_data,
    "adj_mtx" = 999
  ))

}
