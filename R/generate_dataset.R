#' Generate sub-dataset based on latin square.
#'
#' A function which extracts the corresponding subdataset based on a latin
#' square.
#'
#' @param number integer vector; latin square numbers. Length is equal to the
#' number of returned data frames.
#' @param ls_mat matrix; latin square index matrix.
#' @param adjacency matrix; network adjacency matrix.
#' @param adjacency2 matrix; squared network adjacency matrix.
#' @param use.multicore logical; if set to \code{TRUE} (and if package
#' \code{\link{parallel}} is available) the function is run in parallel,
#' default is \code{FALSE}.
#' @param mc.cores count; number of cores to use for parallel execution if
#' \code{use.multicore = FALSE}, default is 2.
#'
#' @return A list of data frames is returned, where each data frame has columns
#' \describe{
#'  \item{y}{binary response.}
#'  \item{edges}{edges change statictic (always 1).}
#'  \item{twostars}{2-star change statictic.}
#'  \item{triangles}{triangle change statictic.}
#' }
#'
#' @seealso \code{\link{ergam}}
#'
#' @examples
#' library("network")
#' data("florentine", package = "ergm")
#'
#' ls_mat <- simpleLS(n = get.network.attribute(flomarriage, "n"))
#' adj_mat <- as.matrix(flomarriage)
#'
#' test_dat <- generate_dataset(number = 1,
#'                              ls_mat = ls_mat,
#'                              adjacency = adj_mat,
#'                              adjacency2 = adj_mat %*% adj_mat)
#' test_dat
#'
generate_dataset <- function(number = 1,
                             ls_mat,
                             adjacency,
                             adjacency2,
                             use.multicore = FALSE,
                             mc.cores = 2)
{
  FUN <- function(number,
                  adjacency,
                  adjacency2) {

    sampling_indices <- which(ls_mat == number, arr.ind = TRUE)
    sampling_indices <-
      sampling_indices[(sampling_indices[, 1] < sampling_indices[, 2]),]

    ## Extract information for current subsample ###############################

    ## Response vector
    y <- adjacency[sampling_indices]

    x_edges <- rep(1, times = length(y))

    ## Compute 2-star change
    k <- 2
    deg_seq <- diag(adjacency2)
    x_2stars <- choose((deg_seq[sampling_indices[, 1]] - 1)^y +
                         (deg_seq[sampling_indices[, 1]])^(1 - y) - 1, k - 1) +
      choose((deg_seq[sampling_indices[, 2]] - 1)^y +
               (deg_seq[sampling_indices[, 2]])^(1 - y) - 1, k - 1)

    ## Extract triangle change
    x_triangles <- adjacency2[sampling_indices]

    X <- data.frame(y = y,
                    edges = x_edges,
                    twostars = x_2stars,
                    triangles = x_triangles)

    rownames(X) <- NULL

    attr(X, "ls_num") <- number

    return(X)
  }

  if(use.multicore & require("parallel")) {
    data_list <- mclapply(number,
                          FUN = FUN,
                          adjacency = adjacency,
                          adjacency2 = adjacency2,
                          mc.cores = mc.cores)
  } else {
    data_list <- lapply(number,
                        FUN = FUN,
                        adjacency = adjacency,
                        adjacency2 = adjacency2)
  }

  return(data_list)
}
