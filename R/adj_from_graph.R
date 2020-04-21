#' Generate a graph corresponding to the "lattice" graph
#'
#' @param graph Graph, represented as a data frame with two columns ('n1' and
#'     'n2'). Each row represents an edge between two nodes. Values in the
#'     data frame represent the index of the node
#' @return An adjacency matrix
#' @export
adj_from_graph <- function(graph) {

  dim <- max(c(graph$n1,graph$n2))
  adj_mtx <- matrix(0, nrow=dim, ncol=dim)

  for (i in 1:nrow(graph)) {

    n1 <- graph[i,"n1"]
    n2 <- graph[i,"n2"]
    adj_mtx[n1,n2] <- 1
    adj_mtx[n2,n1] <- 1

  }

  return (adj_mtx)

}
