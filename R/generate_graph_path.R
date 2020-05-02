#' Generate a graph corresponding to the "path" graph
#'
#' @description Numbering starts on the left (#1) and goes to the right (#100)
#' @return A data frame representing the edges of a graph
#' @export
generate_graph_path <- function() {

  data <- data.frame(
    n1 = 1:99,
    n2 = 2:100
  )

  return (data)

}
