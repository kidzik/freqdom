#' Converts a multivariate time series or a time-domain filter into a time-domain object, enabling further manipulation using
#' time domain methods such as \code{\link{fourier.transform}} or \code{\link{filter.process}}.
#'
#' \code{timedom} is a technical class of objects which simplifies manipulation of lagged multivariate linear filters and multivariate time series. It
#' enables representation of a series of filters (a set of matrices) in a single object, which can be used in
#' \code{\link{fourier.transform}}, \code{\link{filter.process}} and other functions.
#' For generality and simplicity of representation, we assume that a multivariate time series \eqn{X_t}
#' with elements in \eqn{\mathbf{R}^p} and \eqn{t \in \{1,2,...T\}} is also a set of matrices of size \eqn{p \times 1}.
#' With such representation, functions such as \code{\link{fourier.transform}} are invariant to the type
#' of the object, i.e. they behave the same for time-domain filters and for multivariate time series.
#'
#' For a given set of operators \eqn{\{B_k : k \in lags\}}, such that \eqn{X_k \in \mathbf{R}^{p_1\times p_2}}
#' for \code{lags} \eqn{ \subset \mathbf{Z}}, the function \code{timedom} creates an object which we refer to as a \code{timedom} object.
#' 
#' If the data is given as an array \code{X} of size \eqn{L \times p_1 \times p_2} and a vector of \code{lags} is
#' provided, then the element \code{X[i,,]} corresponds tothe lag \code{lags[i]} (see example). If \code{lags} are not given
#' we assume \code{lags} be \eqn{\{1,2,...,L\}}. 
#' 
#' If the data is given as a matrix \eqn{X} of size \eqn{T \times p_1} representing a time series,
#' then we assume it's of dimensions \eqn{T \times p_1 \times 1} and set \code{lags} to \eqn{\{1,2,...,T\}}. 
#'
#' @title Converts a multivariate time series or a time-domain filter into a time-domain object
#' 
#' @param X a set of \eqn{L} operators represented as an array of size \eqn{L \times p_1 \times p_2} or
#' a matrix of size \eqn{T \times p} representing a multivariate time series
#' @param lags a vector of length \eqn{L} of lags on which operators are defined. \code{lags[i]} corresponds to the operator \code{X[i,,]}
#' @return A time domain object represented as a list with elements
#' * \code{$operators} - defined by \eqn{X},
#' * \code{$lags} - \code{lags} if they are provided or, otherwise, \eqn{\{1,2,...,L\}} where \eqn{L} is the number of observations in \eqn{X},
#' 
#' dependnig on the type of the input.
#' @export
#' 
#' @examples
#' d = 2
#' 
#' A = array(0,c(d,d,2))
#' A[1,,] = 2 * diag(d:1)/d
#' A[2,,] = 1.5 * diag(d:1)/d
#' OP = timedom(A,c(-2,1))
#' print(OP)
timedom = function (X,lags=NULL){
  res = list()

  if (is.vector(X)){
    if (is.null(lags))
      lags = 1:length(X)
    res$operators = array(0,c(length(X),1,1))
    res$operators[,1,1] = X
  }
  else if (is.matrix(X)){
    if (is.null(lags))
      lags = 1:dim(X)[1]

    res$operators = array(0,c(dim(X)[1],dim(X)[2],1))
    res$operators[,,1] = X
  }
  else if (is.array(X) && length(dim(X)) == 3){
    if (is.null(lags))
      lags = 1:dim(X)[3] - 1
    res$operators = X
  }
  else{
    stop("X must be a matrix or an array")
  }
  
  res$lags = lags
  class(res) = 'timedom'
  res
}