#' @include build_basis.R
NULL

#' Visualize basis functions of ergam.
#'
#' Function to plot the used basis functions for an \code{\link{ergam}} fit.
#'
#' @param m \code{\link{ergam}} object.
#' @param xlim numeric vector; plotting range on x-axis. If not specified,
#' plotting range is computed automatically.
#'
#' @seealso \code{\link{ergam}}
#'
#' @export
#'
visualize_basis <- function(m,
                            xlim = NULL) {

  if(is.null(xlim)) {
    xlim <- c(0, max(m$max_twostar_obs, m$max_triangle_obs))
  }

  x_seq <- seq(from = xlim[1], to = xlim[2])

  if (!is.matrix(m$parameters_values)) {
    parameters_values <- as.matrix(m$parameters_values, ncol = 1)
  } else {
    parameters_values <- m$parameters_values
  }

  distfun <- get(paste("p", m$disttype, sep = ""),
                 mode = "function", envir = parent.frame())

  plot(x = x_seq, y = FUN_dist(x = x_seq,
                               par_value = parameters_values[1, ],
                               distfun = distfun,
                               parameters_names = m$parameters_names),
       xlim = xlim, type = "l", ylim = c(0, 1),
       main = "Used basis functions for smooth ergam fit",
       xlab = "x", ylab = "Value of basis function")

  for (i in 2:nrow(parameters_values)) {
    lines(x = x_seq, y = FUN_dist(x = x_seq,
                                  par_value = parameters_values[i,],
                                  distfun = distfun,
                                  parameters_names = m$parameters_names))
  }
}
