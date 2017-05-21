#' Given a time-domain filter \eqn{\{A_k : k \in S\}} (\code{\link{timedom}}) for some lags \eqn{S}, \code{timedom.norms}
#' returns \eqn{\{\|A_k\| : k \in S\}} for a given norm \eqn{\| \cdot \|}.
#' 
#' @title Compute operator norms of elements of a filter
#' @param A a time-domain filter \code{\link{timedom}}, i.e. a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
#' for some set \eqn{S \subset \mathbf{Z}} of lags.
#' @param type the matrix norm to be used as in \link{base::norm}
#' @export
#' @return A list with elements
#' * \code{$norms} - a vector of the length \eqn{S} with the \eqn{k}-th covariate equal \eqn{\|A_{lags(k)}\|},
#' * \code{$lags} - \code{lags} corresponding to the computed norms.
#' @examples
#' d = 2
#' 
#' A = array(0,c(d,d,2))
#' A[1,,] = 2 * diag(d:1)/d
#' A[2,,] = 1.5 * diag(d:1)/d
#' OP = timedom(A,c(-2,1))
#' timedom.norms(OP)
timedom.norms = function(A, type='2'){
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
    v = c(v,mynorm(matrix(A$operators[i,,],D[2],D[3])))
  R$norms = v
  R
}

