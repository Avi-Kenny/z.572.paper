#' Generate an adjacency matrix corresponding to the "United States" graph
#'
#' @return An adjacency matrix
#' @export
adj_from_gis_us <- function() {

  # Adapted from https://www.counterpointstat.com/uploads/1/1/9/3/
  #              119383887/areal.r

  usa_state <- map(database="state", fill=TRUE, plot=FALSE)
  state_id <- sapply(strsplit(usa_state$names, ":"), function(x) x[1])
  usa_poly <- map2SpatialPolygons(usa_state, IDs=state_id)
  usa_nb <- poly2nb(usa_poly)
  adj_mtx <- nb2mat(usa_nb, style="B")
  adj_mtx[(adj_mtx>0)] <- 1
  adj_mtx <- adj_mtx[c(1:7,9:49),c(1:7,9:49)]
  attr(adj_mtx, "dimnames") <- NULL

  return (adj_mtx)

}
