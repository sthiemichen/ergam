#' @include ergam_methods.R
NULL

#' Plot ergam results.
#'
#' Function to plot resulting estimates for \code{\link{ergam}}.
#'
#' @param m \code{\link{ergam}} object.
#' @param what character string; specifies what results should be plotted.
#' Default is \code{"converged"}. See details.
#' @param xlim numeric vector; plotting range on x-axis. If not specified,
#' plotting range is computed automatically.
#' @param plot_rug logical; should the observed values be plotted? Default is
#' \code{FALSE}. \cr
#' Advise: Do not set to \code{TRUE} for big networks.
#' @param show_mean logical; should a mean curve be plotted (in blue)? Default
#' is \code{TRUE}.
#' @param show_median logical; should a median curve be plotted (in orange)?
#' Default is \code{TRUE}.
#' @param alpha_shade numeric; value between 0 and 1 specifying the alpha
#' transparency for the plotted curve estimates.
#'
#' @details
#' The default option is to plot only converged estimates (including estimates
#' where smooth effects have been set to zero). If \code{what = "all"} is used,
#' all estimated effects are plotted (for non-converged estimates, the last
#' value from the iterative algorithm is used).
#' Option \code{what = "non-converged"} may be useful for debugging or
#' diagnostics if a lot of subsets encounter convergence problems. \cr
#' If a median curve is dislayed (\code{show_median = TRUE}) the function
#' \code{\link{fbplot}} from package \code{\link{fda}}
#' (Ramsay et. al., 2014) is used with argument \code{method = "Both"}.
#' The computation here relies on the fast algorithm developed by
#' Sun et. al. (2012).
#'
#' @seealso \code{\link{ergam}}
#'
#' @references
#' Ramsay, J. O., Wickham, H., Graves, S. and Hooker, G. (2014).
#' "fda: Functional Data Analysis". R package version 2.4.4.
#' \url{https://CRAN.R-project.org/package=fda} \cr\cr
#' Sun, Y., Genton, M. G. and Nychka, D. W. (2012),
#' "Exact fast computation of band depth for large functional datasets:
#' How quickly can one million curves be ranked?" Stat, 1(1), 68-74.
#'
#' @examples
#' \dontrun{
#' data("facebook")
#'
#' mod1 <- ergam(nw = facebook, ls_num = 1:2)
#'
#' plot(mod1)
#' }
#'
#' @export
#'
plot.ergam <- function(m,
                       what = c("converged", "all", "nonconverged"),
                       xlim = NULL,
                       plot_rug = FALSE,
                       show_mean = TRUE,
                       show_median = TRUE,
                       alpha_shade = 0.5) {

  .pardefault <- par(no.readonly = TRUE)

  what <- match.arg(what)

  if(plot_rug) {
    data_temp <- whole_net_to_data(m$nw)
  }

  res <- generate_full_mat(m = m,
                           what = what,
                           med = show_median)

  if(m$mod == "full") {
    par(mfrow = c(3, 1))
  } else {
    par(mfrow = c(2, 1))
  }

  P <- boxplot(res$f_all_mat[, 1],
               ylab = expression(theta[edges]))
  title(main = paste("Estimates for intercept (",
                     nrow(res$f_all_mat), " estimates)", sep = ""))

  segments(1:length(P$n) - 0.1976, P$stats[3, ], 1:length(P$n) + 0.1976,
           P$stats[3, ], lwd = 5, col = "White") # hide marginal median

  if(show_mean) {
    abline(h = mean(res$f_all_mat[, 1]), lwd = 2, col = "red")
  }

  if(show_median) {
    abline(h = res$f_all_mat[res$median_index, 1], lwd = 2, col = "orange")
  }

  if(res$par_2) {

    if(is.null(xlim)) {

      pred_min <- 0
      pred_max <- m$max_twostar_obs

    } else {

      pred_min <- xlim[1]
      pred_max <- xlim[2]

    }

    plot(res$pred_seq_twostar, res$f_all_mat[1, res$u2_f_index], type = "l",
         xlab = expression(Delta[twostar]),
         ylab = expression(m(Delta[twostar])), cex = 1.2,
         xlim = c(pred_min, pred_max),
         ylim = range(res$f_all_mat[, res$u2_f_index]),
         main = paste("Estimates for smooth twostar effect (",
                      nrow(res$alpha_mat), " estimates)", sep = ""),
         col = rgb(0.8, 0.8, 0.8, alpha = alpha_shade))

    if(plot_rug) {
      rug(data_temp$twostars)
    }

    for (i in 2:nrow(res$f_all_mat)) {
      lines(res$pred_seq_twostar, res$f_all_mat[i, res$u2_f_index],
            col = rgb(0.8, 0.8, 0.8, alpha = alpha_shade))
    }

    if(show_mean) {
      lines(res$pred_seq_twostar, apply(res$f_all_mat[, res$u2_f_index],
                                        MARGIN = 2, mean),
            lwd = 2, col = "blue")
    }

    if(show_median) {
      lines(res$pred_seq_twostar,
            res$f_all_mat[res$median_index, res$u2_f_index],
            lwd = 2, col = "orange")
    }

  }

  if(res$par_3) {

    if(is.null(xlim)) {

      pred_min <- 0
      pred_max <- m$max_triangle_obs

    } else {

      pred_min <- xlim[1]
      pred_max <- xlim[2]

    }

    plot(res$pred_seq_triangle, res$f_all_mat[1, res$u3_f_index], type = "l",
         xlab = expression(Delta[triangle]),
         ylab = expression(m(Delta[triangle])), cex = 1.2,
         xlim = c(pred_min, pred_max),
         ylim = range(res$f_all_mat[, res$u3_f_index]),
         main = paste("Estimates for smooth triangle effect (",
                      nrow(res$alpha_mat), " estimates)", sep = ""),
         col = rgb(0.8, 0.8, 0.8, alpha = alpha_shade))

    if(plot_rug) {
      rug(data_temp$triangles)
    }

    for (i in 2:nrow(res$f_all_mat)) {
      lines(res$pred_seq_triangle, res$f_all_mat[i, res$u3_f_index],
            col = rgb(0.8, 0.8, 0.8, alpha = alpha_shade))
    }

    if(show_mean) {
      lines(res$pred_seq_triangle, apply(res$f_all_mat[, res$u3_f_index],
                                         MARGIN = 2, mean),
            lwd = 2, col = "blue")
    }

    if(show_median) {
      lines(res$pred_seq_triangle,
            res$f_all_mat[res$median_index, res$u3_f_index],
            lwd = 2, col = "orange")
    }

  }

  par(.pardefault)

}

#' @export
#'
plot.ergam_residuals <- function(resid,
                                 ...) {

  plot(resid$residuals,
       main = paste("Mean Pearson residuals per node based on",
                    resid$estimate_type,  "ergam model"), ...)

}
