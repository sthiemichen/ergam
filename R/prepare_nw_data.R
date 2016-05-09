#' Prepare network data for ergam.
#'
#' A function to extract the needed information for the \code{\link{ergam}}
#' function from the whole network data set.
#'
#' @param adjacency matrix; the network adjacency matrix.
#'
#' @details Note: If the number of nodes in the network is odd, the generated
#' latin square has dimension (no. of nodes - 1).
#'
#' @return A list with elements
#' \describe{
#'  \item{nnodes}{numeric; number of nodes in the network.}
#'  \item{max_2star_change}{numeric; maximal value for the 2-star
#'  change.}
#'  \item{max_triangle_change}{numeric; maximal value for the triangle
#'  change.}
#'  \item{adjacency}{matrix; network adjacency matrix.}
#'  \item{adjacency2}{matrix; squared network adjacency matrix.}
#'  \item{ls_mat}{matrix; corresponding latin square.}
#' }
#'
#' @seealso \code{\link{ergam}}
#'
#' @examples
#' library("network")
#' data("florentine", package = "ergm")
#'
#' test_nw_dat <- prepare_nw_data(as.matrix(flomarriage))
#' test_nw_dat
#'
prepare_nw_data <- function(adjacency) {
  nnodes <- ncol(adjacency)

  max_twostar_change <- 2*nnodes - 2
  max_triangle_change <- nnodes - 2

  if (nnodes %% 2 == 0) { # For Latin square, n needs to be even
    ls_mat <- simpleLS(n = nnodes)
  } else {
    ls_mat <- simpleLS(n = nnodes - 1)
  }

  res <- list(nnodes = nnodes,
              max_twostar_change = max_twostar_change,
              max_triangle_change = max_triangle_change,
              adjacency = adjacency,
              adjacency2 = crossprod(adjacency),
              ls_mat = ls_mat)
  return(res)
}
