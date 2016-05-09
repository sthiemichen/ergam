#' @include prepare_nw_data.R
NULL
#' @include build_basis.R
NULL
#' @include generate_dataset.R
NULL
#' @include fit_mod.R
NULL
#' @include ergam_methods.R
NULL

#' Exponential Random Graph Additive Model.
#'
#' Fitting an Exponential Random Graph Additive Model (ERGAM) with smooth
#' components.
#'
#' @param nw \code{\link{network}} object or adjacency \code{\link{matrix}}.
#' Network has to be undirected.
#' @param mod character string; type of model to fit. Should be one of \cr
#' \code{"full"} for model with 2-star and triangle effect, \cr
#' \code{"twostar"} for model with only 2-star effect, or \cr
#' \code{"triangle"} for model with only triangle effect. \cr
#' Default is \code{"full"}.
#' @param sig integer (vector), contains the direction of the smooth effects:
#' \code{-1} for decreasing, \code{1} for increasing.
#' Two entries needed for full model (first entry for 2-star, second for
#' triangle effect). \cr
#' If not specified (default), the effect direction is determined using a GLM
#' fit for each data subset.
#' @param disttype character string; type of basis functions to use for smooth
#' effect(s), see details. \cr
#' Default is \code{"exp"} for the exponential distribution.
#' @param parameters_names character (vector); containing the names for the
#' distribution parameters. \cr
#' Default is \code{"rate"} for the exponential distribution.
#' @param parameters_values numeric vector or matrix; containing a sequence of
#' values for the distribution parameters. The no. of columns must be equal to
#' the no. of parameters. See details.
#' @param nodes numeric; positions where monotonicity constraint is checked,
#' see details.
#' @param base_num integer; number of basis functions to use, see details.
#' @param ls_num integer (vector); specifies which latin square based subsets
#' should be used for model fitting.\cr
#' The default value of \code{0} indicates all available subsets.
#' @param iter_max integer; maximum number of iterations, default is 20.
#' @param y_min integer; minimum number of observations with response 1. If a
#' subset contains less, fitting is skipped for this subset.
#' @param exact numeric; convergence criterion, default is \code{1e-12}.
#' @param use.multicore logical; should the algorithm run in parallel mode, if
#' package \link{parallel} is available? Default is \code{FALSE}.
#' @param mc.cores integer; no. of cores to use for parallel execution, default
#' is 2.
#' @param check logical; prints some model diagnostics using function
#' \code{\link{check.ergam}}. Default is \code{TRUE}.
#' @param stab_fac numeric; factor for avoiding numerical issues in
#' \code{\link{quadprog}}, default is \code{1e-9}.
#' @param zero_tolerance numeric; tolerance for setting effect to zero in case
#' of numeric difficulties. Default is \code{0.005}.
#'
#' @details
#' The function fits an Exponential Random Graph Additive Model (ERGAM)
#' with smooth components to subsets of the whole network. The subsets are
#' chosen based on a symmetric Latin square. \cr
#'
#' The argument \code{disttype} specifies the distribution type for the basis
#' functions of the smooth effect(s). Anything can be used where the cumulative
#' distribution function is available with function name \code{pdisttype} and
#' the argument for the quantiles named \code{q}, and the corresponding density
#' function as \code{disttype} with argument name \code{x} for the quantile
#' vector. Most standard distributions in \R satisfy this condition. \cr
#' If no \code{parameter_values} are sepcified and \code{"exp"} is used as
#' \code{disttype} (default), 20 basis functions are used for argument
#' \code{base_num}. In addition, for the exponential distribution the
#' monotonicity constraint is checked at points where neighbouring basis
#' functions cross each other (\code{nodes} argument). \cr
#'
#' The actual fitting for \code{\link{ergam}} with smooth
#' twostar and/or triangle effect is done using the \code{\link{quadprog}}
#' routines (Turlach and Weingessel, 2013) to incorporate the
#' monotonicity constraints (direction as specified using \code{sig} argument or
#' based on automatic GLM fits per data subset).\cr
#' A smooth effect (twostar or triangle) can be set to zero for a data subset
#' during fitting (happens when the corresponding penalty parameter tends to
#' infinity, or if numeric issues occur in a later than the first iteration
#' round and the effect does not exceed the value specified in
#' \code{zero_tolerance}).
#'
#' @return an object of class \code{ergam}, for which
#'   \code{\link[=plot.ergam]{plot}},
#'   \code{\link[=check.ergam]{check}},
#'   \code{\link[=residuals.ergam]{residuals}}, and
#'   \code{\link[=gof.ergam]{gof}} methods are defined. Objects of class
#'   \code{ergam} are lists with entries
#'   \describe{
#'      \item{\code{results}}{a list containing the results for each of the data
#'      subset used, with entries \cr
#'      \code{alpha} for the estimated coefficients, \cr
#'      \code{lambda} for the estimated penalization parameter(s), \cr
#'      \code{status} containing information on the algorithm convergence, \cr
#'      \code{iter} for the number of iterations used, and \cr
#'      \code{sig} effect direction(s) used.}
#'      \item{\code{nw}}{the supplied \code{\link{network}} object.}
#'      \item{\code{mod}}{type of the fitted model, see argument \code{mod}.}
#'      \item{\code{sig}}{the direction used for the smooth effect(s) (\code{-1}
#'      implies decreasing, \code{1} increasing).}
#'      \item{\code{disttype}}{ype of basis functions used for smooth
#'      effect(s).}
#'      \item{\code{parameters_names}}{names for the distribution parameters for
#'      the basis functions.}
#'      \item{\code{parameters_values}}{values for the distribution distribution
#'      parameters of the basis functions.}
#'      \item{\code{nodes}}{positions where monotonicity constraints have been
#'      checked during fitting.}
#'      \item{\code{base_num}}{no. of basis functions used.}
#'      \item{\code{ls_num}}{information which subsets have been used for model
#'      fitting. \code{0} implies all.}
#'      \item{\code{iter_max}}{maximum no. of iterations.}
#'      \item{\code{y_min}}{minimum no. of observations with response 1 per
#'      subset.}
#'      \item{\code{exact}}{convergence criterion.}
#'      \item{\code{zero_tolerance}}{criterion for setting effect to zero.}
#'      \item{\code{max_twostar_obs}}{max. observed twostar change in data.}
#'      \item{\code{max_triangle_obs}}{max. observed triangle change in data.}
#'      }
#'
#' @references
#' Thiemichen, S. and Kauermann, G. (2016). Stable exponential random graph
#' models with non-parametric components for large dense networks. \cr
#' arXiv preprint arXiv:1604.04732. \url{http://arxiv.org/abs/1604.04732} \cr
#' \cr
#' Turlach, B. A. and Weingessel, A. (2013).
#' "quadprog: Functions to solve Quadratic Programming Problems". \cr
#' R package version 1.5-5. \url{https://CRAN.R-project.org/package=quadprog}
#'
#' @examples
#' \dontrun{
#' # Example akes some time to be executed
#' data("facebook")
#'
#' mod1 <- ergam(nw = facebook, ls_num = 1:10)
#'
#' plot(mod1)
#'
#' plot(residuals(mod1))
#' }
#'
#' @import network
#'
#' @export ergam
#'
ergam <- function(nw,
                  mod = c("full", "twostar", "triangle"),
                  sig = NULL,
                  disttype = "exp",
                  parameters_names = "rate",
                  parameters_values = NULL,
                  nodes = NULL,
                  base_num = NULL,
                  ls_num = 0,
                  iter_max = 20,
                  y_min = 10,
                  exact = 1e-12,
                  use.multicore = FALSE,
                  mc.cores = 2,
                  check = TRUE,
                  stab_fac = 1e-9,
                  zero_tolerance = 0.005) {

  if(!(class(nw) == "network")) {
    if(class(nw) == "matrix") {
      nw <- network(nw, directed = FALSE)
    } else {
      stop("Please use either a network object or an adjacency matrix!")
    }
  }

  if(is.directed(nw)) {
    stop("Network has to be undirected.")
  }

  mod <- match.arg(mod)

  if(is.null(sig)) { # automatic effect direction specification ################
    if(mod == "full") {
      sig_val <- c(1, 1)
    } else {
      sig_val <- 1
    }
  } else { # handle user input for effect direction ############################
    if(mod == "full") {
      if(length(sig) == 1) {
        sig_val <- c(sig, sig)
      } else {
        sig_val <- sig
      }
    } else {
      if (length(sig) == 1) {
        sig_val <- sig
      } else {
        stop("Only one effect direction can be specified in twostar and triangle model!")
      }
    }
  }


  nw_data <- prepare_nw_data(adjacency = as.matrix(nw))


  if ((length(ls_num) == 1) && (ls_num == 0)) {
    if (nw_data$nnodes %% 2 == 0) { # For Latin square, n needs to be even
      ls_num <- seq(from = 1, to = (nw_data$nnodes - 1), by = 1)
    } else {
      ls_num <- seq(from = 1, to = (nw_data$nnodes - 2), by = 1)
    }
  }


  data_list <- generate_dataset(number = ls_num,
                                adjacency = nw_data$adjacency,
                                adjacency2 = nw_data$adjacency2,
                                ls_mat = nw_data$ls_mat,
                                use.multicore = use.multicore,
                                mc.cores = mc.cores)

  max_obs <- apply(matrix(unlist(lapply(data_list,
                                        FUN = function(u) {
                                          return(c(max(u$twostars),
                                                   max(u$triangles)))
                                        })),
                          ncol = 2, byrow = TRUE),
                   MARGIN = 2, max)
  max_twostar_obs <- max_obs[1]
  max_triangle_obs <- max_obs[2]

  if(is.null(parameters_values)) {
    if(disttype == "exp") {
      if(is.null(base_num)) { base_num <- 20 }
      max_change_obs  <- max(unlist(data_list))
      parameters_values <- 1 / exp(seq(from = log(max_change_obs), to = 0,
                                       length.out = base_num))
    } else {
      stop("Please specify values for the distribution parameters.")
    }
  } else if(is.vector(parameters_values)) {
    parameters_values <- sort(parameters_values)
    base_num <- length(parameters_values)
  } else {
    base_num <- nrow(parameters_values)
  }

  if(is.null(nodes)) {
    if(disttype == "exp") {
      nodes <- (log(parameters_values[-length(parameters_values)]) -
                  log(parameters_values[-1])) /
        (parameters_values[-length(parameters_values)] - parameters_values[-1])
    } else {
      stop("Nodes need to be specified manually, if the exponential distribution is not employed.")
    }
  }

  if (mod == "full") {
    X_current2 <- build_basis(x = nodes, disttype = disttype,
                              parameters_values = parameters_values,
                              parameters_names = parameters_names,
                              dens = TRUE)
    X_current3 <- build_basis(x = nodes, disttype = disttype,
                              parameters_values = parameters_values,
                              parameters_names = parameters_names,
                              dens = TRUE)
    A <- rbind(cbind(rep(0, times = nrow(X_current2)),
                     sig_val[1] * X_current2,
                     matrix(0, nrow = nrow(X_current2),
                            ncol = ncol(X_current3))),
               cbind(rep(0, times = nrow(X_current3)),
                     matrix(0, nrow = nrow(X_current3),
                            ncol = ncol(X_current2)),
                     sig_val[2] * X_current3))

    D2 <- diag(c(0, rep(c(1, 0), each = (ncol(A) - 1) / 2)))
    D3 <- diag(c(0, rep(c(0, 1), each = (ncol(A) - 1) / 2)))
    D <- NULL

  } else {

    X_current <- build_basis(x = nodes, disttype = disttype,
                             parameters_values = parameters_values,
                             parameters_names = parameters_names,
                             dens = TRUE)
    A <- cbind(rep.int(0, times = nrow(X_current)), sig_val * X_current)

    D <- diag(c(0, rep(1, each = ncol(A) - 1)))
    D2 <- NULL
    D3 <- NULL

  }

  if(use.multicore & require("parallel")) {
    res <- mclapply(data_list,
                    FUN = fit_mod,
                    A = A,
                    D = D,
                    D2 = D2,
                    D3 = D3,
                    disttype = disttype,
                    parameters_values = parameters_values,
                    parameters_names = parameters_names,
                    nodes = nodes,
                    mod = mod,
                    sig = sig,
                    iter_max = iter_max,
                    y_min = y_min,
                    exact = exact,
                    stab_fac = stab_fac,
                    zero_tolerance = zero_tolerance,
                    mc.cores = mc.cores)
  } else {
    res <- lapply(data_list,
                  FUN = fit_mod,
                  A = A,
                  D = D,
                  D2 = D2,
                  D3 = D3,
                  disttype = disttype,
                  parameters_values = parameters_values,
                  parameters_names = parameters_names,
                  mod = mod,
                  sig = sig,
                  iter_max = iter_max,
                  y_min = y_min,
                  exact = exact,
                  stab_fac = stab_fac,
                  # subset_id = TRUE, # for debugging
                  zero_tolerance = zero_tolerance)

  }

  results <- list(results = res,
                  nw = nw,
                  mod = mod,
                  sig = sig,
                  disttype = disttype,
                  parameters_names = parameters_names,
                  parameters_values = parameters_values,
                  nodes = nodes,
                  base_num = base_num,
                  ls_num = ls_num,
                  iter_max = iter_max,
                  y_min = y_min,
                  exact = exact,
                  zero_tolerance = zero_tolerance,
                  max_twostar_obs = max_twostar_obs,
                  max_triangle_obs = max_triangle_obs)

  class(results) <- "ergam"

  if (check) {
    check(results, print_only = TRUE)
  }

  return(results) ##############################################################

}
