#' Setup basis functions.
#'
#' Setup basis functions based on selected distribution type.
#'
#' @param x numeric vector.
#' @param disttype character string; see details.
#' Default is \code{"exp"} for the exponential distribution.
#' @param parameters_values numeric vector or matrix; containing a sequence of
#' values for the distribution parameters. The no. of columns must be equal to
#' the no. of parameters. \cr
#' Default is a sequence from 1 to 5.
#' @param parameters_names character vector; containing the names for the
#' distribution parameters. \cr
#' Default is \code{"rate"} for the exponential distribution.
#' @param dens logical; if \code{TRUE} the density values are returned instead
#' of the cumulative distribution function.
#' Default is \code{FALSE}.
#'
#' @details The argument \code{disttype} specifies the distribution type.
#' Anything can be used where the cumulative distribution function is available
#' with function name \code{pdisttype} and the argument for the
#' quantiles named \code{q}, and the corresponding density function as
#' \code{ddisttype} with argument name \code{x} for the quantile vector.
#' Most standard distributions in \R satisfy this condition.
#'
#' @return A matrix is returned containing the values of the cumulative
#' distribution function (or density). Each column represents one parameter
#' configuration, that is, one row from the \code{parameters_values} argument.
#'
#' @examples
#' x <- 1:10
#' base_mat <- build_basis(x = x)
#' base_mat
#'
build_basis <- function(x,
                        disttype = "exp",
                        parameters_values = seq(from = 1, to = 5, by = 1),
                        parameters_names = c("rate"),
                        dens = FALSE) {
  if (is.character(disttype)) {
    densfun <- get(paste("d", disttype, sep = ""),
                   mode = "function", envir = parent.frame())
    distfun <- get(paste("p", disttype, sep = ""),
                   mode = "function", envir = parent.frame())
  } else {
    stop("argument disttype needs to be a chracter string!")
  }

  if (!is.matrix(parameters_values)) {
    parameters_values <- as.matrix(parameters_values, ncol = 1)
  }

  if(dens) {
    return(apply(parameters_values, MARGIN = 1, FUN = FUN_dens,
                 x = x, densfun = densfun, parameters_names))
  } else {
    return(apply(parameters_values, MARGIN = 1, FUN = FUN_dist,
                 x = x, distfun = distfun, parameters_names))
  }
}


FUN_dist <- function(par_value, x, distfun, parameters_names) {
  args <- append(list(x), as.list(par_value))
  names(args) <- c("q", parameters_names)

  distfun_values <- do.call(distfun, args = args)

  ## basis needs to start at 0 for x = 0
  args[[1]] <- 0
  distfun0 <- do.call(distfun, args = args)

  return(distfun_values - distfun0)
}


FUN_dens <- function(par_value, x, densfun, parameters_names) {
  args <- append(list(x), as.list(par_value))
  names(args) <- c("x", parameters_names)

  densfun_values <- do.call(densfun, args = args)

  return(densfun_values)
}
