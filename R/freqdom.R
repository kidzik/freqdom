#' Converts an array of filters into a frequency-domain object enabling further manipulation using
#' frequency-domain methods such as \code{\link{fourier.inverse}}.
#' Typically, these objects are created through a Fourier transform (\code{\link{fourier.transform}})
#' of time-domain objects (multivariate time series or filters \code{\link{timedom}}). 
#'
#' \code{freqdom} is a technical class of objects which simplifies manipulation of linear operators in frequency domain. It
#' enables representation of a sequence of frequency domain operators (a sequence of matrices) in a single object, which can be used in
#' \code{\link{fourier.inverse}} or any of the \code{freqdom.*} functions.
#'
#' For a given set of operators \eqn{\{X_k : k \in freq \}}, such that \eqn{X_k \in \mathbf{C}^{p_1\times p_2}}
#' for some \code{freq} \eqn{ \subset \mathbf{Z}}, the function \code{freqdom} creates an object which we refer to as a \code{freqdom} object.
#' We assume that frequencies \code{freq} are defined on the interval \eqn{[-\pi,\pi]}.
#' 
#' For an array \code{F} of dimensions \eqn{L \times p_1 \times p_2}, each element \code{F[i,,]} corresponds to freq \code{freq[i]} for \eqn{1 \leq i \leq L}.
#' 
#' @title Converts an array of filters into a frequency-domain object
#' 
#' @param F a set of \eqn{L} operators represented as an array of dimensions \eqn{L \times p_1 \times p_2}.
#' @param freq a vector of length \eqn{L} of frequencies on which the operators are defined. \code{freq[i]} corresponds to the operator \code{F[i,,]}.
#' @return A frequency domain operator represented as a list with elements
#' * \code{$operators} = \eqn{F},
#' * \code{$freq} = \code{freq}.
#' @seealso \code{\link{fourier.transform}}
#' @export
freqdom = function (F,freq)
{
  if (!is.array(F) || length(dim(F)) != 3)
    stop("F must be an array of evaluations")
  
  res = list()
  res$operators = F
  res$freq = freq
  class(res) = 'freqdom'
  res
}