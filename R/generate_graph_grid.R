#' Generate a graph corresponding to the "grid" graph
#'
#' @description Numbering starts on the top left (#1) and goes to the top right
#'     (#10), then continues on the next row (#11), etc.
#' @return A data frame representing the edges of a graph
#' @export
generate_graph_grid <- function() {

  data <- data.frame(
    "n1" = integer(),
    "n2" = integer()
  )

  for (i in 1:100) {

    # Add left-right neighbors
    if (i%%10 != 0) {
      data[nrow(data)+1,] <- c(i,i+1)
    }

    # Add up-down neighbors
    if (i %in% c(1:90)) {
      data[nrow(data)+1,] <- c(i,i+10)
    }

  }

  return (data)

}
