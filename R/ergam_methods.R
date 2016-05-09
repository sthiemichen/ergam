#' @include build_basis.R
NULL

generate_full_mat <- function(m,
                              what,
                              med) {

  check_results <- check(m, print_only = FALSE, no_print = TRUE)

  all_alpha <- lapply(m$results, FUN = function (u) {u$alpha})

  alpha_temp <- switch(what,
                       converged = all_alpha[(check_results$which_converged |
                                                check_results$which_effect_zero)],
                       all = all_alpha[!(check_results$which_no_fit)],
                       nonconverged =
                         all_alpha[check_results$which_non_convergence])

  alpha_mat <- matrix(unlist(alpha_temp),
                      nrow = length(alpha_temp), byrow = TRUE)

  where_u2 <- logical(length = ncol(alpha_mat))
  where_u3 <- logical(length = ncol(alpha_mat))

  if(m$mod == "full") {

    par_2 <- TRUE
    par_3 <- TRUE

    where_u2[2:((ncol(alpha_mat) - 1)/2 + 1)] <- TRUE
    where_u3[((ncol(alpha_mat) - 1)/2 + 2):ncol(alpha_mat)] <- TRUE

  } else if (m$mod == "triangle")  {

    par_2 <- FALSE
    par_3 <- TRUE

    where_u3[-1] <- TRUE

  } else {

    par_2 <- TRUE
    par_3 <- FALSE

    where_u2[-1] <- TRUE

  }

  pred_min <- 0

  f_all_mat <- matrix(alpha_mat[, 1], ncol = 1) # intercept values

  u2_f_index <- NULL
  u3_f_index <- NULL

  if(par_2) { # twostar

    pred_seq_twostar <- seq(from = pred_min,
                            to = m$max_twostar_obs, by = 1)

    X_pred <- build_basis(x = pred_seq_twostar, disttype = m$disttype,
                          parameters_values = m$parameters_values,
                          parameters_names = m$parameters_names)

    pred_values <- t(apply(alpha_mat[, where_u2], MARGIN = 1,
                           FUN = function(u) {X_pred %*% u}))

    u2_f_index <- seq(from = (ncol(f_all_mat) + 1),
                      to = (ncol(pred_values) + 1))

    f_all_mat <- cbind(f_all_mat, pred_values)
  }

  if(par_3) { # triangle

    pred_seq_triangle <- seq(from = pred_min,
                             to = m$max_triangle_obs, by = 1) # triangle

    X_pred <- build_basis(x = pred_seq_triangle, disttype = m$disttype,
                          parameters_values = m$parameters_values,
                          parameters_names = m$parameters_names)

    pred_values <- t(apply(alpha_mat[, where_u3], MARGIN = 1,
                           FUN = function(u) {X_pred %*% u}))  # triangle

    u3_f_index <- seq(from = (ncol(f_all_mat) + 1),
                      to = (ncol(f_all_mat) + ncol(pred_values)))

    f_all_mat <- cbind(f_all_mat, pred_values)
  }

  if (med) {
    median_index <- get_median(f_all_mat)
  } else {
    median_index <- NULL
  }

  return(list(f_all_mat = f_all_mat,
              alpha_mat = alpha_mat,
              median_index = median_index,
              u2_f_index = u2_f_index,
              u3_f_index = u3_f_index,
              where_u2 = where_u2,
              where_u3 = where_u3,
              par_2 = par_2,
              par_3 = par_3,
              pred_seq_twostar = pred_seq_twostar,
              pred_seq_triangle = pred_seq_triangle))

}

#' Model diagnostics for ergam.
#'
#' Function to summarize algorithm convergence for \code{\link{ergam}} and
#' print some diagnostics.
#'
#' @param m \code{\link{ergam}} object.
#'
#' @export
check <- function(m, ...) UseMethod("check", m)


#' @param print_only logical; diagnostic information are only printed.
#' Default is \code{TRUE}.
#'
#' @return Nothing is returned, unless option \code{print_only = FALSE} is
#' specified. In this case a list is returned with elements
#'   \describe{
#'      \item{\code{which_converged}}{logical vector; indicates for which data
#'      subsets the algorithm did converge.}
#'      \item{\code{which_no_fit}}{logical vector; indicates for which data
#'      subsets there was no fit (due to too few observations with \eqn{y = 1}
#'      in the subset).}
#'      \item{\code{which_max_iter}}{logical vector; indicates for which data
#'      subsets the algorithm did reach the max. no. of iterations.}
#'      \item{\code{which_non_convergence}}{logical vector; indicates for which
#'      data subsets the algorithm did not converge.}
#'      \item{\code{which_effect_zero}}{logical vector; indicates for which data
#'      subsets at least one smooth effect was set to zero.}
#'      }
#'
#' @seealso \code{\link{ergam}}
#'
#' @examples
#' \dontrun{
#' data("facebook")
#'
#' ## By default ergam function includes check functionality
#' mod1 <- ergam(nw = facebook, ls_num = 1, check = TRUE)
#'
#' check(mod1)
#' }
#'
#' @export
#'
#' @rdname check
check.ergam <- function(m,
                        print_only = TRUE,
                        no_print = FALSE) {

  all_status <- lapply(m$results, FUN = function (u) {u$status})

  ## find status == "converged"
  which_converged <- all_status == "converged"

  ## find status == "no fit"
  which_no_fit <- all_status == "no fit"

  ## find status == "maximum number of iterations"
  which_max_iter <- all_status == "maximum number of iterations"
  ## find status == "non-convergence"
  which_non_convergence <- all_status == "non-convergence"

  ## find status == "effect zero"
  which_effect_zero <- all_status == "effect zero"

  if(!no_print) {
    cat("Model diagnostics for ergam: \n")
    cat("Total no. of subsets used:", length(m$results), "\n")
    cat("No. of subsets with convergence: ", sum(which_converged), "\n")
    cat("No. of subsets with no fit (less than", m$y_min,
        "times y = 1 in subset): ", sum(which_no_fit), "\n")
    cat("No. of subsets where max. no. iterations was reached (change max.iter for fit, if needed):",
        sum(which_max_iter), "\n")
    cat("No. of subsets with other reason for non-convergence:",
        sum(which_non_convergence), "\n")

    if (m$mod == "full") {
      if (any(which_effect_zero)) {
        effect_zero_lambda_mat <- matrix(unlist(lapply(m$results[which_effect_zero],
                                                       FUN = function (u) {u$lambda})),
                                         ncol = 2, byrow = TRUE)
        two_NA <- is.na(effect_zero_lambda_mat[, 1])
        tri_NA <- is.na(effect_zero_lambda_mat[, 2])
        two <- sum(two_NA & !tri_NA)
        tri <- sum(!two_NA & tri_NA)
        both <- sum(two_NA & tri_NA)
        cat("No. of subsets where only smooth twostar effect was set to zero:",
            two, "\n")
        cat("No. of subsets where only smooth triangle effect was set to zero:",
            tri, "\n")
        cat("No. of subsets where both smooth effects were set to zero:",
            both, "\n")
      } else {
        cat("No. of subsets where smooth effect was set to zero: 0 \n")
      }
    } else {
      cat("No. of subsets where smooth effect was set to zero:",
          sum(which_effect_zero), "\n")
    }
  }

  if(!print_only) {
    return(list(which_converged = which_converged,
                which_no_fit = which_no_fit,
                which_max_iter = which_max_iter,
                which_non_convergence = which_non_convergence,
                which_effect_zero = which_effect_zero))
  }
}

get_median <- function(f_mat) {
  if (!requireNamespace("fda", quietly = TRUE)) {
    stop("Please install package fda: install.packages('fda')")
  }

  return(suppressWarnings(fda::fbplot(t(f_mat),
                                      method = "Both", plot = FALSE)$medcurve))
}

#' Median estimate for ergam.
#'
#' Function to get median estimate from an \code{\link{ergam}} model.
#'
#' @param m \code{\link{ergam}} object.
#' @param what character string; specifies what should be used for overall
#' estimate calculation. Default is \code{"converged"}. See details.
#'
#' @details
#' The default option is to consider only converged estimates (including
#' estimates where smooth effects have been set to zero) for calculation of
#' the overall median estimate.
#' If \code{what = "all"} is used, all estimated effects are used
#' (for non-converged estimates, the last value from the iterative
#' algorithm is used). \cr
#' The function \code{\link{fbplot}} from package \code{\link{fda}}
#' (Ramsay et. al., 2014) is used with argument \code{method = "Both"}.
#' The computation here relies on the fast algorithm developed by
#' Sun et. al. (2012). \cr
#'
#' @return a named list with entries
#'   \describe{
#'      \item{\code{edges}}{numeric; the model intercept.}
#'      \item{\code{m_twostars}}{numeric vector; containing the values of the
#'      smooth functional curve for the twostar change starting at a change of
#'      0.}
#'      \item{\code{m_triangles}}{numeric vector; containing the values of the
#'      smooth functional curve for the triangle change starting at a change of
#'      0.}
#'   }
#'
#' @seealso \code{\link{ergam}}.
#'
#' @import network
#'
#' @export
#'
median.ergam <- function(m,
                         what = c("converged", "all", "nonconverged")) {

  what <- match.arg(what)

  res <- generate_full_mat(m = m,
                           what = what,
                           med = TRUE)

  nnodes <- get.network.attribute(m$nw, "n")

  seq_twostars <- seq(from = 0, to = (2 * nnodes - 4))
  seq_triangles <- seq(from = 0, to = (nnodes - 2))

  alpha_est <- res$alpha_mat[res$median_index, ]

  if(m$mod == "triangle") {
    alpha_est <- c(alpha_est[1], rep(0, times = m$base_num), alpha_est[-1])
  } else if(m$mod == "twostar") {
    alpha_est <- c(alpha_est, rep(0, times = m$base_num))
  }

  X_twostars <- build_basis(x = seq_twostars,
                            disttype = m$disttype,
                            parameters_values = m$parameters_values,
                            parameters_names = m$parameters_names)

  X_triangles <- build_basis(x = seq_triangles,
                             disttype = m$disttype,
                             parameters_values = m$parameters_values,
                             parameters_names = m$parameters_names)
  results <- list(edges = alpha_est[1],
                  m_twostars = X_twostars %*% alpha_est[2:(m$base_num + 1)],
                  m_triangles = X_triangles %*% alpha_est[-(1:(m$base_num + 1))])

  return(results)
}

#' Mean estimate for ergam.
#'
#' Function to get mean estimate from an \code{\link{ergam}} model.
#'
#' @param m \code{\link{ergam}} object.
#' @param what character string; specifies what should be used for overall
#' estimate calculation. Default is \code{"converged"}. See details.
#'
#' @details
#' The default option is to consider only converged estimates (including
#' estimates where smooth effects have been set to zero) for calculation of
#' the overall mean estimate.
#' If \code{what = "all"} is used, all estimated effects are used
#' (for non-converged estimates, the last value from the iterative
#' algorithm is used).
#'
#' @return a named list with entries
#'   \describe{
#'      \item{\code{edges}}{numeric; the model intercept.}
#'      \item{\code{m_twostars}}{numeric vector; containing the values of the
#'      smooth functional curve for the twostar change starting at a change of
#'      0.}
#'      \item{\code{m_triangles}}{numeric vector; containing the values of the
#'      smooth functional curve for the triangle change starting at a change of
#'      0.}
#'   }
#'
#' @seealso \code{\link{ergam}}.
#'
#' @import network
#'
#' @export
#'
mean.ergam <- function(m,
                       what = c("converged", "all", "nonconverged")) {

  what <- match.arg(what)

  res <- generate_full_mat(m = m,
                           what = what,
                           med = TRUE)

  nnodes <- get.network.attribute(m$nw, "n")

  seq_twostars <- seq(from = 0, to = (2 * nnodes - 4))
  seq_triangles <- seq(from = 0, to = (nnodes - 2))

  alpha_est <- apply(res$alpha_mat, MARGIN = 2, mean)

  if(m$mod == "triangle") {
    alpha_est <- c(alpha_est[1], rep(0, times = m$base_num), alpha_est[-1])
  } else if(m$mod == "twostar") {
    alpha_est <- c(alpha_est, rep(0, times = m$base_num))
  }

  X_twostars <- build_basis(x = seq_twostars,
                            disttype = m$disttype,
                            parameters_values = m$parameters_values,
                            parameters_names = m$parameters_names)

  X_triangles <- build_basis(x = seq_triangles,
                             disttype = m$disttype,
                             parameters_values = m$parameters_values,
                             parameters_names = m$parameters_names)
  results <- list(edges = alpha_est[1],
                  m_twostars = X_twostars %*% alpha_est[2:(m$base_num + 1)],
                  m_triangles = X_triangles %*% alpha_est[-(1:(m$base_num + 1))])

  return(results)
}

#' Residuals for ergam.
#'
#' Function to compute mean Pearson residuals for each node based on
#' \code{\link{ergam}} model.
#'
#' @param m \code{\link{ergam}} object.
#' @param what character string; specifies what should be used for overall
#' estimate calculation. Default is \code{"converged"}. See details.
#' @param estimate_type character string; specifies which estimate should be
#' used for residual computation, \code{"mean"} (default) or \code{"median"}
#' curve.
#'
#' @details
#' An overall estimate is used to compute Pearson residuals for all in the
#' network adjacency matrix. For each node in the network (i.e for each row in
#' the adjacency) the mean value of these residuals is computed. \cr
#' The default option is to employ only converged estimates (including estimates
#' where smooth effects have been set to zero) for calculation of
#' the overall estimate (mean or median, see argument \code{estimate_type}).
#' If \code{what = "all"} is used, all estimated effects are used
#' (for non-converged estimates, the last value from the iterative
#' algorithm is used). \cr
#' If a median curve is used as overall estimate
#' (\code{estimate_type = "median"}) the function
#' \code{\link{fbplot}} from package \code{\link{fda}}
#' (Ramsay et. al., 2014) is used with argument \code{method = "Both"}.
#' The computation here relies on the fast algorithm developed by
#' Sun et. al. (2012). \cr
#'
#' @return an object of class \code{ergam_residuals}, for which a
#'   \code{\link[=plot.ergam_residuals]{plot}} method is defined.
#'   Objects of class
#'   \code{ergam_residuals} are lists with entries
#'   \describe{
#'      \item{\code{residuals}}{numeric vector; containing the mean Pearson
#'      residuals for each node.}
#'      \item{\code{what}}{character string; which results were used for overall
#'      estimate calculation.}
#'      \item{\code{estimate_type}}{character string; which estimate was used
#'      for goodness-of-fit diagnostics (mean or median curve).}
#'      }
#'
#' @seealso \code{\link{ergam}}.
#'
#' @import network
#'
#' @export
#'
residuals.ergam <- function(m,
                            what = c("converged", "all", "nonconverged"),
                            estimate_type = c("mean", "median")) {
  what <- match.arg(what)
  estimate_type <- match.arg(estimate_type)

  res <- generate_full_mat(m = m,
                           what = what,
                           med = (estimate_type == "median"))

  nw_full <- whole_net_to_data(m$nw)

  nnodes <- get.network.attribute(m$nw, "n")
  adjacency <- as.matrix(m$nw)

  alpha_est <- switch(estimate_type,
                      mean = apply(res$alpha_mat, MARGIN = 2, mean),
                      median = res$alpha_mat[res$median_index, ])

  if(m$mod == "triangle") {
    alpha_est <- c(alpha_est[1], rep(0, times = m$base_num), alpha_est[-1])
  } else if(m$mod == "twostar") {
    alpha_est <- c(alpha_est, rep(0, times = m$base_num))
  }

  X_twostars <- build_basis(x = nw_full$twostars,
                            disttype = m$disttype,
                            parameters_values = m$parameters_values,
                            parameters_names = m$parameters_names)

  X_triangles <- build_basis(x = nw_full$triangles,
                             disttype = m$disttype,
                             parameters_values = m$parameters_values,
                             parameters_names = m$parameters_names)

  pi_mat <- matrix(NA, ncol = ncol(adjacency), nrow = nrow(adjacency))
  indices <- which(upper.tri(adjacency, diag = FALSE), arr.ind = TRUE)

  pi_mat[indices] <- exp(cbind(nw_full$edges, X_twostars, X_triangles) %*%
                           alpha_est) /
    (1 + exp(cbind(nw_full$edges, X_twostars, X_triangles) %*% alpha_est))

  pi_mat[lower.tri(pi_mat)] <- t(pi_mat)[lower.tri(pi_mat)]

  pearson_res_mat <- (adjacency - pi_mat) / sqrt(pi_mat * (1 - pi_mat))

  pearson_resid_node <- apply(pearson_res_mat, MARGIN = 1, mean, na.rm = TRUE)

  result <- list(residuals = pearson_resid_node,
                 what = what,
                 estimate_type = estimate_type)

  class(result) <- "ergam_residuals"

  return(result)

}

