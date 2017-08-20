#' This function applies a linear filter to some vector time series.
#'
#' Let \eqn{[X_1,\ldots, X_T]^\prime} be a \eqn{T\times d} matrix corresponding to a vector series \eqn{X_1,\ldots,X_T}. This time series is transformed to the series \eqn{Y_1,\ldots, Y_T}, where
#' \deqn{ Y_t=\sum_{k=-q}^p A_k X_{t-k},\quad t\in\{p+1,\ldots, T-q\}.}
#' The index \eqn{k} of \eqn{A_k} is determined by the lags defined for the time domain object.
#' When index \eqn{t-k} falls outside the domain \eqn{\{1,\ldots, T\}} we set \eqn{X_t=\frac{1}{T}\sum_{k=1}^T X_k}.
#'
#' @title Convolute (filter) a multivariate time series using a time-domain filter
#'
#' @param X vector time series given in matrix form. Each row corresponds to a timepoint.
#' @param A an object of class \code{\link{timedom}}.
#' @return A matrix. Row \eqn{t} corresponds to \eqn{Y_t}.
#' @seealso \code{\link{timedom}}
#' @export
#' @import matrixcalc
#' @keywords time.domain
#' @describeIn filter.process Multivariate convolution (filter) in the time domain
filter.process = function(X, A){
  if (!is.timedom(A))
    stop("A must be timedom object")
  if (!is.matrix(X))
    stop("X must be a matrix")

  ncols = dim(A$operators)[1]
  n = dim(X)[1]
  Y = matrix(0,n,ncols)

  # Filters the process
  for (component in 1:ncols)
  {
    # Compute inner products of X with all lags
    IP = t(X %*% A$operators[component,,])

    # Shift the inner products by lags
    # so that we can just sum them up
    for (i in 1:length(A$lags))
      IP[i,] = colshift(IP[i,], A$lags[i])

    # sum them up!
    Y[,component] = colSums(IP)
  }
  Y
}

# shifts a column by a given lag and imputes mean value of the column for
# values from outside
colshift = function(col,lag){
  n = length(col)
  mcol = mean(col)
  if (abs(lag) >= n)
    rep(mcol,n)
  else if (lag<0)
    c(col[-(1:(-lag))],rep(mcol,(-lag)))
  else if (lag>0)
    c(rep(mcol,lag),col[(lag+1):n - lag])
  else
    col
}

# Given a multivariate time series \eqn{X_t} and a time-domain filter \eqn{\{A_k : k \in S\}} for some lags \eqn{S}, \code{filter.process} convolutes
# the time series with the filter. This is a convenience operator equivalent to \code{\link{filter.process}}. See \code{\link{filter.process}} for details.
#
# @title Convolution of a process X with an operator A.
# @param X multivariate time series represented as a \eqn{T \times p_1} matrix,
#  where \eqn{T} is the number of observations and \eqn{p_1} is the number of covariates.
# @param A time-domain filter \code{\link{timedom}}, i.e. a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
# for some set \eqn{S \subset \mathbf{Z}} of lags.
# @return Multivariate time series \eqn{Y_t} represented as a matrix \eqn{T \times p_2}.
# @seealso \code{\link{filter.process}}, \code{\link{timedom}}
#' @describeIn filter.process Convenience operator for \code{filter.process} function
#' @keywords time.domain
#' @export
`%c%` <- function(X, A){
  filter.process(X, A)
}
