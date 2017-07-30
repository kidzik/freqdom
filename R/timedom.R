#@exportClass timedom
#setClass("timedom", representation(operators = "array", lags = "vector"))

#' Creates an object of class \code{timedom}. This object corresponds to a multivariate linear filter.
#'
#' This class is used to describe a linear filter, i.e. a sequence of matrices, each of which correspond to a certain lag. Filters can, for example, be used to transform a sequence \eqn{(X_t)} into a new sequence \eqn{(Y_t)} by defining
#' \deqn{
#'  Y_t=\sum_k A_kX_{t-k}.
#' }
#' See \code{filter.process()}.
#'  Formally we consider a collection \eqn{[A_1,\ldots,A_K]} of complex-valued matrices \eqn{A_k}, all of which have the same dimension \eqn{d_1\times d_2}. Moreover, we consider lags \eqn{\ell_1<\ell_2<\cdots<\ell_K}. The object this function creates corresponds to the mapping \eqn{f: \mathrm{lags}\to \mathbf{R}^{d_1\times d_2}}, where \eqn{\ell_k\mapsto A_k}.
#'
#' @title Defines a linear filter
#'
#' @param A a vector, matrix or array. If array, the elements \eqn{A[,,k], 1\leq k\leq K}, are real valued \eqn{(d_1\times d_2)} matrices (all of same dimension). If A is a matrix, the \eqn{k}-th row is treated as \eqn{A[,,k]}. Same for the \eqn{k}-th element of a vector. These matrices, vectors or scalars define a linear filter.
#' @param lags a vector of increasing integers. It corresponds to the time lags of the filter.
#' @return Returns an object of class \code{timedom}. An object of class  \code{timedom} is a list containing the following components:
#' * \code{operators} \eqn{\quad} returns the array \code{A} as given in the argument.
#' * \code{lags} \eqn{\quad} returns the vector \code{lags} as given in the argument.
#' @seealso \code{\link{freqdom}}, \code{\link{is.timedom}}
#' @export
#' @keywords classes
#'
#' @examples
#' # In this example we apply the difference operator: Delta X_t= X_t-X_{t-1} to a time series
#' X = rar(20)
#' OP = array(0,c(2,2,2))
#' OP[,,1] = diag(2)
#' OP[,,2] = -diag(2)
#' A = timedom(OP, lags = c(0,1))
#' filter.process(X, A)
timedom = function (A,lags)
{
  res = list()

  if (is.vector(A)){
#    if (is.null(lags))
#      lags = 1:length(A)
	if(length(A) != length(lags))
	stop("length of A must match number of lags")
    res$operators = array(0,c(1,1,length(A)))
    res$operators[1,1,] = A
  }
  else if (is.matrix(A)){
#    if (is.null(lags))
#      lags = 1:dim(A)[1]
	if(length(lags) != dim(A)[1])
	stop("number of rows of A must match number of lags")
    res$operators = array(0,c(1,dim(A)[2],dim(A)[1]))
    res$operators[1,,] = t(A)
  }
  else if (is.array(A) && length(dim(A)) == 3)
  {
#    if (is.null(lags))
#      lags = 1:dim(A)[3] - 1
	if(length(lags) != dim(A)[3])
	stop("number of matrices of A must match number of lags")
    res$operators = A
  }
  else{
    stop("A must be a vector, a matrix or an array")
  }
  if(length(lags)==1)
  dimnames(res$operators)[[3]]<-list(paste("lag", lags))
  dimnames(res$operators)[[3]]<-paste("lag", lags)
  res$lags = lags
  class(res) = 'timedom'
  res
}

