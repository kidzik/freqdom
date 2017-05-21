#' Given a multivariate time series \eqn{X_t} and a time-domain filter \eqn{\{A_k : k \in S\}} defined on lags \eqn{S}, \code{filter.process} convolutes
#' the time series with the filter.
#' 
#' Let \eqn{X_t} be a multivariate time series and let \eqn{\{A_k : k \in S\}} be a time domain filter defined on lags \eqn{S}
#' For a sample of size \eqn{T} and for \eqn{t \in \{1,...,T\}}, the functions \code{filter.process} computes
#' \deqn{Y_t = \sum_{k \in \mathbf{Z}} A_k X_{t-k}}
#' assuming \eqn{X_t = 0} for \eqn{t \notin \{1,...,T\}} and \eqn{A_k = 0} for \eqn{k \notin S}.
#' 
#' @title Convolute (filter) a multivariate time series using a time-domain filter
#' 
#' @param X a multivariate time series represented with a \eqn{T \times p_1} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_1} is the number of covariates.
#' @param A a time-domain filter \code{\link{timedom}}, i.e. a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
#' defined on a set of lags \eqn{S \subset \mathbf{Z}}. 
#' @return A multivariate time series \eqn{Y_t} represented as a matrix \eqn{T \times p_2}.
#' @seealso \code{\link{timedom}}
#' @export
filter.process = function(X, A){
  if (!is.timedom(A))
    stop("A must be timedom object")
  if (!is.matrix(X))
    stop("X must be a multivariate time series (a matrix of observations)")
  
  n = dim(X)[1]
  nbasis = dim(A$operators)[2]
  Y = matrix(0,n,nbasis)
  
  for (component in 1:nbasis)
  {
    IP = A$operators[,component,] %*% t(X)
    for (i in 1:length(A$lags)){
      lag = A$lags[i]
      IP[i,] = colshift(IP[i,],lag)
    }
    Y[,component] = colSums(IP)
    
  }
  Y
}

colshift = function(col,lag){
  n = length(col)
  if (abs(lag) >= n)
    rep(0,n)
  else if (lag<0)
    c(col[-(1:(-lag))],rep(0,(-lag)))
  else if (lag>0)
    c(rep(0,lag),col[(lag+1):n - lag])
  else
    col
}


#' Given a multivariate time series \eqn{X_t} and a time-domain filter \eqn{\{A_k : k \in S\}} for some lags \eqn{S}, \code{filter.process} convolutes
#' the time series with the filter. This is a convenience operator equivalent to \code{\link{filter.process}}. See \code{\link{filter.process}} for details.
#'  
#' @title Convolution of a process X with an operator A.
#' @param X multivariate time series represented as a \eqn{T \times p_1} matrix,
#'  where \eqn{T} is the number of observations and \eqn{p_1} is the number of covariates.
#' @param A time-domain filter \code{\link{timedom}}, i.e. a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
#' for some set \eqn{S \subset \mathbf{Z}} of lags. 
#' @return Multivariate time series \eqn{Y_t} represented as a matrix \eqn{T \times p_2}.
#' @seealso \code{\link{filter.process}}, \code{\link{timedom}}
#' @export
`%c%` <- function(A, X){
  filter.process(X, A)
}
