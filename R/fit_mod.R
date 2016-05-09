#' @include build_basis.R
NULL

#' Fit ergam model.
#'
#' Function to do the actual fitting of an \code{\link{ergam}} model.
#'
#' @param data_sub data frame; generated using \code{\link{generate_dataset}}.
#' @param A matrix; defines the additional monotonicity constraints taken into
#' account by \code{\link{quadprog}} routines.
#' @param D matrix; penalty matrix in \code{"twostar"} or \code{"triangle"}
#' model.
#' @param D2 matrix; penalty matrix for twostar effect in \code{"full"} model.
#' @param D3 matrix; penalty matrix for triangle effect in \code{"full"} model.
#' @param disttype character string; type of basis functions to use for smooth
#' effects, see details. \cr
#' Default is \code{"exp"} for the exponential distribution.
#' @param parameters_names character (vector); containing the names for the
#' distribution parameters. \cr
#' Default is \code{"rate"} for the exponential distribution.
#' @param parameters_values numeric vector or matrix; containing a sequence of
#' values for the distribution parameters. The no. of columns must be equal to
#' the no. of parameters. See details.
#' @param nodes numeric; positions where monotonicity constraint is checked,
#' see details.
#' @param base_num integer; number of basis functions to use. See details.
#' @param mod character string; type of model to fit. Should be one of \cr
#' \code{"full"} for model with twostar and triangle effect, \cr
#' \code{"twostar"} for model with only twostar effect, or \cr
#' \code{"triangle"} for model with only triangle effect. \cr
#' Default is \code{"full"}.
#' @param sig integer (vector), contains the direction of the effects: \code{-1}
#' for decreasing, \code{1} for increasing. Two entries needed for full model
#' (first entry for twostar, second for triangle effect). \cr
#' If not specified (default), the effect direction is determined using a GLM
#' fit.
#' @param iter_max integer; maximum number of iterations, default is 20.
#' @param y_min integer; minimum number of observations with response 1.
#' @param exact numeric; convergence criterion, default is \code{1e-12}.
#' @param verbose logical; printing of current values? Default is \code{FALSE}.
#' @param subset_id logical; printing of current latin square number? Default is
#' \code{FALSE}.
#' @param stab_fac numeric; factor for avoiding numerical issues in
#' \code{\link{quadprog}}, default is \code{1e-9}.
#' @param zero_tolerance numeric; tolerance for setting effect to zero in case
#' of numeric difficulties. Default is \code{0.005}.
#'
#' @details
#' The argument \code{disttype} specifies the distribution type for the basis
#' functions of the smooth effect(s). Anything can be used where the cumulative
#' distribution function is available with function name \code{pdisttype} and
#' the argument for the quantiles named \code{q}, and the corresponding density
#' function as \code{ddisttype} with argument name \code{x} for the quantile
#' vector. Most standard distributions in \R satisfy this condition. \cr
#' If no \code{parameter_values} are specified and \code{"exp"} is used as
#' \code{disttype} (default), 20 basis functions are used for argument
#' \code{base_num}. In addition, for the exponential distribution the
#' \code{nodes} argument is set so that the monotonicity constraint is checked
#' at points where neighbouring basis functions cross each other.\cr
#' If numeric issues occur in a later than the first iteration round, one smooth
#' effect (or both) is set to zero if it does not exceed the value specified in
#' \code{zero_tolerance}.
#'
#' @return A list is returned, with elements
#' \describe{
#'  \item{alpha}{coefficient values.}
#'  \item{lambda}{penalty parameter(s). Contains two entries (\code{lambda2} and
#'  \code{lambda3}) for a full model.}
#'  \item{status}{information on algorithm status.}
#'  \item{iter}{number of iterations executed.}
#' }
#' See returned objects of \code{\link{fit_double}} for a full model, and
#' \code{\link{fit_single}} otherwise.
#'
#' @seealso \code{\link{fit_double}}, \code{\link{fit_single}}
#'
fit_mod <- function(data_sub,
                    A,
                    D = NULL,
                    D2 = NULL,
                    D3 = NULL,
                    disttype = "exp",
                    parameters_names = "rate",
                    parameters_values = NULL,
                    nodes = NULL,
                    base_num = NULL,
                    mod = "full",
                    sig = NULL,
                    iter_max = 20,
                    y_min = 10,
                    exact = 1e-12,
                    verbose = FALSE,
                    subset_id = FALSE,
                    stab_fac = 1e-9,
                    zero_tolerance = 0.005) {

  if (subset_id) { # helpful for debugging
    cat(attr(data_sub, "ls_num"), "\n")
  }

  y <- data_sub$y

  if(sum(y) < y_min) {
    return(return(list(alpha = NA,
                       lambda = NA,
                       status = "no fit",
                       iter = 0)))
  }

  if(is.null(sig)) {
    if (mod == "full") {
      helpmod_coef <- glm(y ~ twostars + triangles, data = data_sub,
                          family = binomial(link = "logit"))$coef[-1]
      sig <- sign(helpmod_coef)
      if(any(sig == -1)) {
        A <- A %*% diag(c(0, rep(sig, each = (ncol(A) - 1) / 2)))
      }
    } else {
      if(mod == "twostar") {
        mod_formula <- 'y ~ twostars'
      } else {
        mod_formula <- 'y ~ triangles'
      }

      helpmod_coef <- glm(mod_formula, data = data_sub,
                          family = binomial(link = "logit"))$coef[-1]
      if(sign(helpmod_coef) == -1) {
        A <- (-1) * A
      }
    }
  }

  if(is.null(parameters_values) & (disttype == "exp")) {
    if(is.null(base_num)) { base_num <- 20 }
    parameters_values <- 1 / exp(seq(from = log(max(data_sub)), to = 0,
                                     length.out = base_num))
  }

  if(mod == "full") {
    par_2 <- TRUE
    par_3 <- TRUE
  } else if (mod == "triangle")  {
    par_2 <- FALSE
    par_3 <- TRUE
  } else {
    par_2 <- TRUE
    par_3 <- FALSE
  }

  if (par_2) {
    X_twostars <- build_basis(x = data_sub$twostars,
                              disttype = disttype,
                              parameters_values = parameters_values,
                              parameters_names = parameters_names)
  }
  if (par_3) {
    X_triangles <- build_basis(x = data_sub$triangles, disttype = disttype,
                               parameters_values = parameters_values,
                               parameters_names = parameters_names)
  }

  if(mod == "full") {
    C <- cbind(data_sub$edges, X_twostars, X_triangles)
  } else if (par_2) {
    C <- cbind(data_sub$edges, X_twostars)
  } else {
    C <- cbind(data_sub$edges, X_triangles)
  }

  if (mod == "full") {

    res <- fit_double(y = y, D2 = D2, D3 = D3, C = C, A = A,
                      iter_max = iter_max, exact = exact, verbose = verbose,
                      stab_fac = stab_fac, zero_tolerance = zero_tolerance,
                      sig = sig)


  } else {

    res <- fit_single(y = y, D = D, C = C, A = A,
                      iter_max = iter_max, exact = exact, verbose = verbose,
                      stab_fac = stab_fac, zero_tolerance = zero_tolerance,
                      sig = sig)

  }

  return(res)

}

#' Fit simple ergam model.
#'
#' Fit simple \code{\link{ergam}} model with smooth twostar or triangle effect.
#'
#' @param y numeric vector; binary response.
#' @param C matrix; design matrix for the model.
#' @param D matrix; penalty matrix for the model.
#' @param A matrix; defines the additional monotonicity constraints taken into
#' account by \code{\link{quadprog}} routines.
#' @param iter_max integer; maximum number of iterations, default is 20.
#' @param exact numeric; convergence criterion, default is \code{1e-12}.
#' @param verbose logical; printing of current values? Default is \code{FALSE}.
#' @param stab_fac numeric; factor for avoiding numerical issues in
#' \code{\link{quadprog}}, default is \code{1e-9}.
#' @param zero_tolerance numeric; tolerance for setting effect to zero in case
#' of numeric difficulties. Default is \code{0.005}.
#' @param sig integer; giving the effect direction used.
#'
#' @details Function executes actual fit for \code{\link{ergam}} with one smooth
#' effect (twostar or triangle) using the \code{\link{quadprog}} routines
#' (Turlach and Weingessel, 2013) to incorporate the additional monotonicity
#' constraints.\cr
#' Non-convergence of the algorithm or other problems during fitting can be seen
#' by the \code{status} entry. In addition, warnings are produced. The function
#' automatically checks if the maximum number of iterations was reached. \cr
#' If numeric issues occur in a later than the first iteration round, the smooth
#' effect is set to zero if it does not exceed the value specified in
#' \code{zero_tolerance}.
#'
#' @return A list is returned, with elements
#' \describe{
#'  \item{alpha}{coefficient values.}
#'  \item{lambda}{penalty parameter.}
#'  \item{status}{information on algorithm status.}
#'  \item{iter}{number of iterations executed.}
#'  \item{sig}{effect direction used.}
#' }
#'
#' @references
#' Turlach, B. A. and Weingessel, A. (2013).
#' "quadprog: Functions to solve Quadratic Programming Problems".
#' R package version 1.5-5. \url{https://CRAN.R-project.org/package=quadprog}
#'
#' @import quadprog
#'
fit_single <- function (y,
                        C,
                        D,
                        A,
                        iter_max = 20,
                        exact = 1e-12,
                        verbose = FALSE,
                        stab_fac = 1e-9,
                        zero_tolerance = 0.005,
                        sig) {

  it_outer <- 1
  alpha_current <- vector(mode = "numeric", length = ncol(C))

  rel_change_lambda <- 1
  lambda_current <- 10

  while ((it_outer <= iter_max) && (rel_change_lambda > exact)) {

    it_inner <- 1

    rel_change_alpha <- 1

    lambdaD_stab <- (stab_fac * lambda_current) * D

    while ((it_inner <= iter_max) && (rel_change_alpha > exact)) {

      alpha_old <- alpha_current
      pi <- exp(C %*% alpha_current) / (1 + exp(C %*% alpha_current))

      qp_sol <- try(
        solve.QP(
          Dmat = t(C * kronecker(matrix(1, 1, ncol(C)),
                                 matrix((pi * (1 - pi)) * stab_fac))) %*% C  +
            lambdaD_stab,
          dvec = t(C) %*% (as.vector(y - pi) * stab_fac) -
            lambdaD_stab %*% alpha_current,
          Amat = t(A),
          bvec = - A %*% alpha_current),
        silent = TRUE)

      if(class(qp_sol) == "try-error") {
        qp_problem <- TRUE
      } else if (any(is.nan(qp_sol$solution))) {
        qp_problem <- TRUE
      } else {
        qp_problem <- FALSE
      }

      if(qp_problem) {

        if(it_inner == 1) {

          max_effect <- max(abs(C[, -1] %*% alpha_old[-1]))
          effect_low <- FALSE
          if(max_effect < zero_tolerance) { effect_low <- TRUE }

          if((it_outer > 1) && effect_low) {

            alpha_current <- vector(mode = "numeric", length = ncol(C))
            alpha_current[1] <- log(sum(y) / (length(y) - sum(y)))

            return(list(alpha = alpha_current,
                        lambda = NA,
                        status = "effect zero",
                        iter = (it_outer - 1),
                        sig = 0))

          } else { # exit

            warning("Algorithm did not converge for some data subset.")
            return(list(alpha = alpha_old,
                        lambda = lambda_current,
                        status = "non-convergence",
                        iter = (it_outer - 1),
                        sig = sig))

          }

        } else {# try to update lambda with last useful value of alpha

          alpha_current <- alpha_old
          it_inner <- iter_max + 1

        }

      } else {

        alpha_current <- alpha_current + qp_sol$solution

        rel_change_alpha <- sqrt(sum((alpha_old - alpha_current)^2)) /
          sqrt(sum(alpha_old^2))

        if(verbose) {

          cat("Relative change at iteration ", it_inner, ": ", rel_change_alpha,
              "\n\n")

        }

        it_inner <- it_inner + 1

      }

    }

    lambda_old <- lambda_current

    pi <- exp(C %*% alpha_current) / (1 + exp(C %*% alpha_current))

    Dmat0 <- (t(C * kronecker(matrix(1, 1, ncol(C)),
                              matrix(pi * (1 - pi)))) %*% C)[-1, -1] # random part only

    Dmat_prod_diag <- try(rowSums(solve(Dmat0 + (lambda_current * D[-1, -1])) *
                                    t(Dmat0)),                                     # diag(solve(Dmat) %*% Dmat0)
                          silent = TRUE)
    if ((class(Dmat_prod_diag) == "try-error") | any(is.nan(Dmat_prod_diag))) {
      warning("Setting effect to zero for some subset.")

      alpha_current <- vector(mode = "numeric", length = ncol(C))
      alpha_current[1] <- log(sum(y) / (length(y) - sum(y)))

      return(list(alpha = alpha_current,
                  lambda = NA,
                  status = "effect zero",
                  iter = (it_outer - 1),
                  sig = 0))

    }

    lambda_current <- as.numeric(sum(Dmat_prod_diag) /
                                   (t(alpha_current[-1]) %*%
                                      alpha_current[-1]))

    if((1 / lambda_current) <= exact) {
      # variance too small, set effect to zero

      alpha_current <- vector(mode = "numeric", length = ncol(C))
      alpha_current[1] <- log(sum(y) / (length(y) - sum(y)))

      return(list(alpha = alpha_current,
                  lambda = NA,
                  status = "effect zero",
                  iter = (it_outer - 1),
                  sig = 0))

    }

    if(verbose) {
      cat("Current value of lambda: ", lambda_current, "\n\n")
    }

    it_outer <- it_outer + 1

  }

  if (it_outer == iter_max) {

    warning("Maximum number of iterations reached.")
    status_message <- "maximum number of iterations"

  } else {

    status_message <- "converged"

  }

  return(list(alpha = alpha_current,
              lambda = lambda_current,
              status = status_message,
              iter = (it_outer - 1),
              sig = sig))
}


#' Fit full ergam model.
#'
#' Fit full \code{\link{ergam}} model with smooth twostar and triangle effect.
#'
#' @param y numeric vector; binary response.
#' @param C matrix; design matrix for the model.
#' @param D2 matrix; penalty matrix for twostar effect.
#' @param D3 matrix; penalty matrix for triangle effect.
#' @param A matrix; defines the additional monotonicity constraints taken into
#' account by \code{\link{quadprog}} routines.
#' @param iter_max integer; maximum number of iterations, default is 20.
#' @param exact numeric; convergence criterion, default is \code{1e-12}.
#' @param verbose logical; printing of current values? Default is \code{FALSE}.
#' @param stab_fac numeric; factor for avoiding numerical issues in
#' \code{\link{quadprog}}, default is \code{1e-9}.
#' @param zero_tolerance numeric; tolerance for setting effect to zero in case
#' of numeric difficulties. Default is \code{0.005}.
#' @param sig integer vector; giving the effect directions used.
#'
#' @details Function executes actual fit for \code{\link{ergam}} with smooth
#' twostar and triangle effect (full model) using the \code{\link{quadprog}}
#' routines (Turlach and Weingessel, 2013) to incorporate the additional
#' monotonicity constraints.\cr
#' Non-convergence of the algorithm or other problems during fitting can be seen
#' by the \code{status} entry. In addition, warnings are produced. The function
#' automatically checks if the maximum number of iterations was reached. \cr
#' If one of the two smooth effects (twostar or triangle) is set to zero
#' (i.e. corresponding penalty parameter tends to infinity), the function
#' \code{\link{fit_single}} is used accordingly to fit an \code{\link{ergam}}
#' model containing only the remaining smooth effect.\cr
#' If numeric issues occur in a later than the first iteration round, one smooth
#' effect (or both) is set to zero if it does not exceed the value specified in
#' \code{zero_tolerance}.
#'
#' @return A list is returned, with elements
#' \describe{
#'  \item{alpha}{coefficient values.}
#'  \item{lambda}{penalty parameters (two entries: \code{lambda2} and
#'  \code{lambda3}).}
#'  \item{status}{information on algorithm status.}
#'  \item{iter}{number of iterations executed.}
#'  \item{sig}{effect directions used.}
#' }
#'
#' @references
#' Turlach, B. A. and Weingessel, A. (2013).
#' "quadprog: Functions to solve Quadratic Programming Problems".
#' R package version 1.5-5. \url{https://CRAN.R-project.org/package=quadprog}
#'
#' @import quadprog
#'
fit_double <- function (y,
                        C,
                        D2,
                        D3,
                        A,
                        iter_max = 20,
                        exact = 1e-12,
                        verbose = FALSE,
                        stab_fac = 1e-9,
                        zero_tolerance = 0.005,
                        sig) {

  it_outer <- 1
  alpha_current <- vector(mode = "numeric", length = ncol(C))

  rel_change_lambda2 <- 1
  rel_change_lambda3 <- 1

  lambda2_current <- 100
  lambda3_current <- 100

  ## requires same number of basis functions ###################################
  where_u2 <- logical(length = ncol(C))
  where_u3 <- logical(length = ncol(C))
  where_u2[2:((ncol(C)-1)/2 + 1)] <- TRUE
  where_u3[((ncol(C) - 1)/2 + 2):ncol(C)] <- TRUE

  while ((it_outer <= iter_max) &&
         ((rel_change_lambda2 > exact) | (rel_change_lambda3 > exact))) {

    it_inner <- 1

    rel_change_alpha <- 1

    lambdaD_stab <- (lambda2_current * stab_fac) * D2 +
      (lambda3_current * stab_fac) * D3

    while ((it_inner <= iter_max) && (rel_change_alpha > exact)) {

      alpha_old <- alpha_current
      pi <- exp(C %*% alpha_current) / (1 + exp(C %*% alpha_current))
      if(any(is.nan(pi))) {

        sign_vec <- sign(exp(C %*% alpha_current)[is.nan(pi)])
        pi[is.nan(pi)] <- 1^(sign_vec == 1) * 0^(sign_vec == -1)

      }

      qp_sol <- try(
        solve.QP(
          Dmat = t(C * kronecker(matrix(1, 1, ncol(C)),
                                 matrix((pi * (1 - pi)) * stab_fac))) %*% C +
                    lambdaD_stab,
          dvec = t(C) %*% (as.vector(y - pi) * stab_fac) -
            lambdaD_stab %*% alpha_current,
          Amat = t(A),
          bvec = - A %*% alpha_current),
        silent = TRUE)
      if(class(qp_sol) == "try-error") {
        qp_problem <- TRUE
      } else if (any(is.nan(qp_sol$solution))) {
        qp_problem <- TRUE
      } else {
        qp_problem <- FALSE
      }
      if(qp_problem) {

        if(it_inner == 1) {

          max_two_effect <- max(abs(C[, where_u2] %*% alpha_old[where_u2]))
          two_low <- FALSE
          if(max_two_effect < zero_tolerance) { two_low <- TRUE }
          max_tri_effect <- max(abs(C[, where_u3] %*% alpha_old[where_u3]))
          tri_low <- FALSE
          if(max_tri_effect < zero_tolerance) { tri_low <- TRUE }

          if((it_outer > 1) && any(c(two_low, tri_low))) {

            if (two_low && (max_two_effect < max_tri_effect)) {
              # set 2-star effect to zero, if it is the smaller one

              res <- fit_single(y = y,
                                D = D3[!(where_u2), !(where_u2)],
                                C = C[ , !(where_u2)], A = A[ , !(where_u2)],
                                iter_max = iter_max, exact = exact,
                                verbose = verbose,
                                stab_fac = stab_fac, sig = sig[1])
              return(list(alpha = c(res$alpha[1],
                                    rep(0, times = sum(where_u2)),
                                    res$alpha[-1]),
                          lambda = c(lambda2 = NA,
                                     lambda3 = res$lambda),
                          status = "effect zero",
                          iter = (it_outer - 1),
                          sig = c(0, res$sig)))

            } else if(tri_low) { # set triangle effect to zero

              res <- fit_single(y = y,
                                D = D2[!(where_u3), !(where_u3)],
                                C = C[ , !(where_u3)], A = A[ , !(where_u3)],
                                iter_max = iter_max, exact = exact,
                                verbose = verbose,
                                stab_fac = stab_fac, sig = sig[2])
              return(list(alpha = c(res$alpha, rep(0, times = sum(where_u3))),
                          lambda = c(lambda2 = res$lambda,
                                     lambda3 = NA),
                          status = "effect zero",
                          iter = (it_outer - 1),
                          sig = c(res$sig, 0)))

            }

          } else { # exit

            warning("Algorithm did not converge for some data subset.")
            return(list(alpha = alpha_old,
                        lambda = c(lambda2 = lambda2_current,
                                   lambda3 = lambda3_current),
                        status = "non-convergence",
                        iter = (it_outer - 1),
                        sig = sig))
          }

        } else {# try to update lambda with last useful value of alpha

          alpha_current <- alpha_old
          it_inner <- iter_max + 1

        }

      } else {

        alpha_current <- alpha_current + qp_sol$solution

        rel_change_alpha <- sqrt(sum((alpha_old - alpha_current)^2)) /
          sqrt(sum(alpha_old^2))

        if(verbose) {
          cat("Relative change at iteration ", it_inner, ": ", rel_change_alpha,
              "\n\n")
        }

        it_inner <- it_inner + 1

      }
    }

    lambda2_old <- lambda2_current
    lambda3_old <- lambda3_current

    pi <- exp(C %*% alpha_current) / (1 + exp(C %*% alpha_current))

    Dmat0 <- (t(C * kronecker(matrix(1, 1, ncol(C)),
                              matrix(pi * (1 - pi)))) %*%
                C)[as.logical(where_u2 + where_u3),
                   as.logical(where_u2 + where_u3)]

    Dmat_prod_diag <- try(
      rowSums(
        solve(Dmat0 +
                lambdaD_stab[as.logical(where_u2 + where_u3),
                             as.logical(where_u2 + where_u3)] * stab_fac^(-1)) *
          t(Dmat0)), # diag(solve(Dmat) %*% Dmat0)
      silent = TRUE)

    if ((class(Dmat_prod_diag) == "try-error") | any(is.nan(Dmat_prod_diag))) {

      if(lambda2_current > lambda3_current) { # set 2-star effect to zero

        res <- fit_single(y = y,
                          D = D3[!(where_u2), !(where_u2)],
                          C = C[ , !(where_u2)], A = A[ , !(where_u2)],
                          iter_max = iter_max, exact = exact, verbose = verbose,
                          stab_fac = stab_fac, sig = sig[1])
        return(list(alpha = c(res$alpha[1],
                              rep(0, times = sum(where_u2)),
                              res$alpha[-1]),
                    lambda = c(lambda2 = NA,
                               lambda3 = res$lambda),
                    status = "effect zero",
                    iter = (it_outer - 1),
                    sig = c(0, res$sig)))

      } else { # set triangle effect to zero

        res <- fit_single(y = y,
                          D = D2[!(where_u3), !(where_u3)],
                          C = C[ , !(where_u3)], A = A[ , !(where_u3)],
                          iter_max = iter_max, exact = exact, verbose = verbose,
                          stab_fac = stab_fac, sig = sig[2])
        return(list(alpha = c(res$alpha, rep(0, times = sum(where_u3))),
                    lambda = c(lambda2 = res$lambda,
                               lambda3 = NA),
                    status = "effect zero",
                    iter = (it_outer - 1),
                    sig = c(res$sig, 0)))

      }

    }

    if(rel_change_lambda2 > exact) {

      lambda2_current <- as.numeric(sum(Dmat_prod_diag[(where_u2[-1])]) /
                                      (t(alpha_current[where_u2]) %*%
                                         alpha_current[where_u2] )) # trace

      if((1 / lambda2_current) <= exact) {
        # variance too small, set 2-star effect to zero

        res <- fit_single(y = y,
                          D = D3[!(where_u2), !(where_u2)],
                          C = C[ , !(where_u2)], A = A[ , !(where_u2)],
                          iter_max = iter_max, exact = exact, verbose = verbose,
                          stab_fac = stab_fac, sig = sig[1])
        return(list(alpha = c(res$alpha[1],
                              rep(0, times = sum(where_u2)),
                              res$alpha[-1]),
                    lambda = c(lambda2 = NA,
                               lambda3 = res$lambda),
                    status = "effect zero",
                    iter = (it_outer - 1),
                    sig = c(0, res$sig)))

      }

    }

    if(rel_change_lambda3 > exact) {

      lambda3_current <- as.numeric(sum(Dmat_prod_diag[(where_u3[-1])]) /
                                      (t(alpha_current[where_u3]) %*%
                                         alpha_current[where_u3]))

      if((1 / lambda3_current) <= exact) {
        # variance too small, set triangle effect to zero

        res <- fit_single(y = y,
                          D = D2[!(where_u3), !(where_u3)],
                          C = C[ , !(where_u3)], A = A[ , !(where_u3)],
                          iter_max = iter_max, exact = exact, verbose = verbose,
                          stab_fac = stab_fac, sig = sig[2])
        return(list(alpha = c(res$alpha, rep(0, times = sum(where_u3))),
                    lambda = c(lambda2 = res$lambda,
                               lambda3 = NA),
                    status = "effect zero",
                    iter = (it_outer - 1),
                    sig = c(res$sig, 0)))

      }

    }

    if(verbose) {
      cat("Current value of lambda2: ", lambda2_current, "\n\n")
      cat("Current value of lambda3: ", lambda3_current, "\n\n")
    }

    it_outer <- it_outer + 1

  }

  if (it_outer == iter_max) {

    warning("Maximum number of iterations reached.")
    status_message <- "maximum number of iterations"

  } else {

    status_message <- "converged"

  }

  return(list(alpha = alpha_current,
              lambda = c(lambda2 = lambda2_current,
                         lambda3 = lambda3_current),
              status = status_message,
              iter = (it_outer - 1),
              sig = sig))

}
