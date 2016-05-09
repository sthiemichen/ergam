#' Facebook data from Stanford Large Network Dataset Collection.
#'
#' This is the combined data from 10 Facebook ego networks from Stanford Large
#' Network Dataset Collection. These data are stored as a undirected
#' \link{network} object containing the edge information only.
#'
#' @format An object of class \code{\link{network}} with 4039 nodes based on the
#' undirected binary adjacency matrix.
#'
#' @references Leskovec, J. and Krevl, A. (2014) SNAP Datasets: Stanford Large
#' Network Dataset Collection. \url{http://snap.stanford.edu/data} \cr
#' McAuley, J. and Leskovec, J. (2012) Learning to Discover Social Circles in
#' Ego Networks. NIPS.
#'
#' @source \url{https://snap.stanford.edu/data/egonets-Facebook.html}
#'
#' @usage data("facebook")
#'
#' @import network
"facebook"
