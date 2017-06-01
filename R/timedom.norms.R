#' This function determines the norms of the matrices defining some linear filter.
#' 
#' Computes \eqn{\|A_h\|} for \eqn{h} in the set of lags belonging to the object A. When type
#' is \code{2} then \eqn{\|A\|} is the spectral radius of \eqn{A}. When type is \code{F}
#' then \eqn{\|A\|} is the Frobenius norm (or the Hilbert-Schmidt norm, or Schatten 2-norm) of
#' \eqn{A}. 
#' 
#' @title Compute operator norms of elements of a filter
#' @param A an object of class \code{timedom}
#' @param type matrix norm to be used as in \link{base::norm}
#' @export
#' @return A list which contains the following components:
#' * \code{lags} a vector containing the lags of \eqn{A}.
#' * \code{norms} a vector containing the norms of the matrices defining \eqn{A}.
#' @examples
#' d = 2
#' 
#' A = array(0,c(d,d,2))
#' A[1,,] = 2 * diag(d:1)/d
#' A[2,,] = 1.5 * diag(d:1)/d
#' OP = timedom(A,c(-2,1))
#' timedom.norms(OP)
timedom.norms = function(A, type = c("O", "I", "F", "M", "2")){
  if (!is.timedom(A))
    stop ("A must be a time domain filter")
  
  mynorm = norm
  if (type=='2')
    mynorm = norm.spec
  R = list()
  R$lags = A$lags
  v = c()
  D = dim(A$operators)
  for (i in 1:D[1])
    v = c(v,mynorm(matrix(A$operators[,,i],D[2],D[3])))
  R$norms = v
  R
}


