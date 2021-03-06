% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergam_methods.R
\name{median.ergam}
\alias{median.ergam}
\title{Median estimate for ergam.}
\usage{
\method{median}{ergam}(m, what = c("converged", "all", "nonconverged"))
}
\arguments{
\item{m}{\code{\link{ergam}} object.}

\item{what}{character string; specifies what should be used for overall
estimate calculation. Default is \code{"converged"}. See details.}
}
\value{
a named list with entries
  \describe{
     \item{\code{edges}}{numeric; the model intercept.}
     \item{\code{m_twostars}}{numeric vector; containing the values of the
     smooth functional curve for the twostar change starting at a change of
     0.}
     \item{\code{m_triangles}}{numeric vector; containing the values of the
     smooth functional curve for the triangle change starting at a change of
     0.}
  }
}
\description{
Function to get median estimate from an \code{\link{ergam}} model.
}
\details{
The default option is to consider only converged estimates (including
estimates where smooth effects have been set to zero) for calculation of
the overall median estimate.
If \code{what = "all"} is used, all estimated effects are used
(for non-converged estimates, the last value from the iterative
algorithm is used). \cr
The function \code{\link{fbplot}} from package \code{\link{fda}}
(Ramsay et. al., 2014) is used with argument \code{method = "Both"}.
The computation here relies on the fast algorithm developed by
Sun et. al. (2012). \cr
}
\seealso{
\code{\link{ergam}}.
}

