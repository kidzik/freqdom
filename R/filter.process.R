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

  lag.min=min(A$lags)
  lag.max=max(A$lags)
  f.rows=dim(A$operators)[1]
  f.cols=dim(A$operators)[2]
  f.lags=length(A$lags)
  f.span=lag.max-lag.min+1
  n=dim(X)[1]
  p=dim(X)[2]
  Y=matrix(rep(0,n*f.rows),ncol=f.rows)

  # we adapt the filter so that all lags between lag.min and lag.max are covered.
  # filter coefficients of missing lags are set equal to zero
  if(f.span != f.lags){
  Anew = array(0,c(f.rows, f.cols, f.span))
  j=1
  for(i in A$lags-lag.min+1){
  Anew[,,i]=A$operators[,,j]
  j=j+1
  }
  A=timedom(Anew,lags=lag.min:lag.max)
  }

  # we put all filter coefficients in one big matrix
  Filter=c()
  if(dim(A$operators)[1]==1){
  for(j in 1:f.span){Filter=cbind(Filter,matrix(A$operators[,,f.span-j+1],nrow=1))}
  } else
  for(j in 1:f.span){
  	Filter=cbind(Filter,A$operators[,,f.span-j+1])
  }
  # we define the observations before index 1 and after index n as empirical mean of X
  V=X
  if(lag.max>0){
  	V=rbind(t(matrix(rep(colMeans(X),lag.max),nrow=p)),X)
  }
  else{
  V=X[(-lag.max+1):(dim(X)[1]),]
  }
  if(lag.min<0){
  	V=rbind(V,t(matrix(rep(colMeans(X),-lag.min),nrow=p)))
  }
  else{
  	V=V[1:(dim(V)[1]-lag.min),]
  }


  for(i in 1:n){
  	Y[i,]=t(Filter%*%vec(t(V[i:(i+f.span-1),])))
  }
  Y
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
