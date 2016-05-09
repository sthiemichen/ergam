#' Simple latin square construction.
#'
#' A function which constructs a simple latin square with \eqn{n} rows and
#' \eqn{n} columns.
#' The latin square is symmetric and all diagonal entries are 0.
#'
#' @param n numeric; dimension of the matrix, needs to be even.
#'
#' @return a symmetric \eqn{n x n} integer matrix with 0 as diagonal elements.
#'
#' @references
#' Bogomolny, A. (2016). Latin squares. Simple construction.
#' From Interactive Mathematics Miscellany and Puzzles
#' \url{http://www.cut-the-knot.org/arithmetic/latin2.shtml}.
#' [Online; accessed: 04 March 2016].
#'
#' @examples
#' testmat <- simpleLS(n = 10)
#' testmat
#'
#' @export simpleLS
#'
simpleLS <- function(n = 10) {
  if(!((n %% 2) == 0)) {
    stop("n has to be even! (At the moment.)")
  }
  if(n < 3) {
    stop("n has to be at least 3!")
  }

  temp <- matrix(NA, nrow = n, ncol = n)

  # First two rows
  temp[1, ] <- c((n-1), 1:(n - 2), NA)
  temp[2, ] <- c(1:(n - 2), (n-1), NA)

  # Remaining rows
  for (i in 3:(n-1)) {
    temp[i, ] <- c((i-1):(n - 2), (n-1), 1:(i-2), NA)
  }

  # Copy diagonal entries to the last row and column
  temp[n, ] <- diag(temp)
  temp[, n] <- diag(temp)

  # Set diagonal entries to 0
  diag(temp) <- rep(0, times = n)

  return(temp)
}
