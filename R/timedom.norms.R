#' This function determines the norms of the matrices defining some linear filter.
#'
#' Computes \eqn{\|A_h\|} for \eqn{h} in the set of lags belonging to the object A. When type
#' is \code{2} then \eqn{\|A\|} is the spectral radius of \eqn{A}. When type is \code{F}
#' then \eqn{\|A\|} is the Frobenius norm (or the Hilbert-Schmidt norm, or Schatten 2-norm) of
#' \eqn{A}. Same options as for the function \code{norm} as in base package.
#'
#' @title Compute operator norms of elements of a filter
#' @param A an object of class \code{timedom}
#' @param type matrix norm to be used as in \code{norm}
#' @export
#' @return A list which contains the following components:
#' * \code{lags} \eqn{\quad} a vector containing the lags of \code{A}.
#' * \code{norms} \eqn{\quad} a vector containing the norms of the matrices defining \code{A}.
#' @keywords time.domain
#' @examples
#' d = 2
#'
#' A = array(0,c(d,d,2))
#' A[1,,] = 2 * diag(d:1)/d
#' A[2,,] = 1.5 * diag(d:1)/d
#' OP = timedom(A,c(-2,1))
#' timedom.norms(OP)
timedom.norms = function(A, type = "2"){
  if (!is.timedom(A))
    stop ("A must be an object of class timedom")
  if (is.vector(type) || type=="2")
  	mytype = "2"
  if (type=="F" || type=="f")
    mytype = "F"
  if (type=="O"|| type=="o")
    mytype = "O"
  if (type=="I" || type=="i")
    mytype = "I"
  if (type=="M"|| type=="m")
    mytype = "M"
  R = list()
  R$lags = A$lags
  v = c()
  D = dim(A$operators)
  for (i in 1:D[3])
    v = c(v,norm(matrix(A$operators[,,i],D[1],D[2]),type=mytype))
  R$norms = v
  R
}


