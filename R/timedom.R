#' Creates an object of class \code{timedom}. This object corresponds to a linear filter.
#'
#' This class is used to describe a linear filter, i.e. a sequence of matrices, each of which correspond to a certain lag. Filters can, for example, be used to transform a sequence \eqn{(X_t)} into a new sequence \eqn{(Y_t)} by defining
#' \deqn{
#'  Y_t=\sum_k A_kX_{t-k}.
#' }
#' See \code{filter.process()}.
#'  Formally we consider a collection \eqn{[A_1,\ldots,A_K]} of complex-valued matrices \eqn{A_k}, all of which have the same dimension \eqn{d_1\times d_2}. Moreover, we consider lags \eqn{\ell_1<\ell_2<\cdots<\ell_K}. The object this function creates corresponds to the mapping \eqn{f: \mathrm{lags}\to \mathbf{R}^{d_1\times d_2}}, where \eqn{\ell_k\mapsto A_k}. 
#'
#' @title Converts a series of matrices into an object which relates a certain lag to each matrix.
#' 
#' @param A The elements \eqn{A[,,k], 1\leq k\leq K}, are real valued \eqn{(d_1\times d_2)} matrices (all of same dimension). These matrices define a linear filter.
#' @param lags a vector of increasing integers. It corresponds to the time lags of the filter.
#' @return Returns an object of class \code{timedom}. An object of class  \code{timedom} is a list containing the following components:
#' * \code{operators} - operators & returns the array A as given in the argument.
#' * \code{lags} - returns the vector of lags as given in the argument.
#' @seealso \code{\link{freqdom}}, \code{\link{is.timedom}}
#' @export
#' 
#' @examples
#' d = 3
#' 
#' A = array(0,c(d,d,2))
#' A[,,1] = 2 * diag(d:1)/d
#' A[,,2] = 1.5 * diag(d:1)/d
#' OP = timedom(A,c(-2,1))
#' print(OP)
timedom = function (A,lags=1:dim(A)[1])
{
  res = list()

  if (is.vector(A)){
    if (is.null(lags))
      lags = 1:length(A)
    res$operators = array(0,c(length(A),1,1))
    res$operators[,1,1] = A
  }
  else if (is.matrix(A)){
    if (is.null(lags))
      lags = 1:dim(A)[1]

    res$operators = array(0,c(dim(A)[1],dim(A)[2],1))
    res$operators[,,1] = A
  }
  else if (is.array(A) && length(dim(A)) == 3)
  {
    if (is.null(lags))
      lags = 1:dim(A)[3] - 1
    res$operators = A
  }
  else{
    stop("A must be a matrix or an array")
  }
  
  res$lags = lags
  class(res) = 'timedom'
  res
}
