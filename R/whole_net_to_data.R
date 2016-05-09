#' Generate dataset from network.
#'
#' A function which extracts the information for an \code{\link{ergam}} fit
#' from the corresponding network.
#'
#' @param nw \code{\link{network}} object.
#' @param range_only logical; should only the \code{\link{range}} of the
#' variables be returned? Default is \code{FALSE}.
#'
#' @return A data frames is returned, with columns
#' \describe{
#'  \item{y}{binary response.}
#'  \item{edges}{edges change statictic (always 1).}
#'  \item{twostars}{2-star change statictic.}
#'  \item{triangles}{triangle change statictic.}
#' }
#' If argument \code{range_only = TRUE} the columns contain only the
#' \code{\link{range}} of the variables.
#'
#' @seealso \code{\link{ergam}}
#'
#' @import network
#'
whole_net_to_data <- function(nw,
                              range_only = FALSE) {

  adjacency <- as.matrix(nw)
  adjacency2 <- crossprod(adjacency)

  indices <- which(upper.tri(adjacency, diag = FALSE), arr.ind = TRUE)

  y <- as.vector(adjacency[indices])

  x_edges <- rep(1, times = length(y))

  ## Compute 2-star change
  k <- 2
  deg_seq <- diag(adjacency2)
  x_2stars <- choose((deg_seq[indices[, 1]] - 1)^y +
                       (deg_seq[indices[, 1]])^(1 - y) - 1, k - 1) +
    choose((deg_seq[indices[, 2]] - 1)^y +
             (deg_seq[indices[, 2]])^(1 - y) - 1, k - 1)

  ## Extract triangle change
  x_triangles <- adjacency2[indices]

  if(range_only) {
    X <- data.frame(y = range(y),
                    edges = range(x_edges),
                    twostars = range(x_2stars),
                    triangles = range(x_triangles))
  } else {
    X <- data.frame(y = y,
                    edges = x_edges,
                    twostars = x_2stars,
                    triangles = x_triangles)
  }
  rownames(X) <- NULL

  return(X)
}
