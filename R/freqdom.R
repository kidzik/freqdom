#@exportClass freqdom
#setClass("freqdom", representation(operators = "array", freq = "vector"))

#' Creates an object of class  \code{freqdom}. This object corresponds to functional with domain \eqn{[-\pi,\pi]} and some complex vector space as codomain.
#'
#' This class is used to describe a frequency domain functional (like a spectral density matrix, a discrete Fourier transform, an impulse response function, etc.)
#' on selected frequencies. Formally we consider a collection \eqn{[F_1,\ldots,F_K]} of complex-valued matrices \eqn{F_k}, all of which have the same dimension
#' \eqn{d_1\times d_2}. Moreover, we consider frequencies \eqn{\{\omega_1,\ldots, \omega_K\}\subset[-\pi,\pi]}. The object this function creates corresponds
#' to the mapping \eqn{f: \mathrm{freq}\to \mathbf{C}^{d_1\times d_2}}, where \eqn{\omega_k\mapsto F_k}.
#' 
#' Consider, for example, the discrete Fourier transform of a vector time series \eqn{X_1,\ldots, X_T}:. It is defined as
#' \deqn{
#'   D_T(\omega)=\frac{1}{\sqrt{T}}\sum_{t=1}^T X_t e^{-it\omega},\quad \omega\in[-\pi,\pi].
#' }
#' We may choose \eqn{\omega_k=2\pi k/K-\pi} and \eqn{F_k=D_T(\omega_k)}. Then, the object \code{freqdom} creates, is corresponding to the function which associates \eqn{\omega_k} and \eqn{D_T(\omega_k)}.
#'
#' @title Converts an array of filters into a frequency-domain object
#' 
#' @param F an array. The elements \eqn{F[,,k], 1\leq k\leq K}, are complex valued \eqn{(d_1\times d_2)} matrices (all of same dimension).
#' @param freq a vector of dimension \eqn{K} containing frequencies in \eqn{[-\pi,\pi]}.
#' @return Returns an object of class \code{\link{freqdom}}. An object of class  \code{\link{freqdom}} is a list containing the following components:
#' * \code{operators} returns the array F as given in the argument.
#' * \code{freq} returns the vector freq as given in the argument.
#' @seealso \code{\link{fourier.transform}}
#' @examples
#' i = complex(imaginary=1)
#' OP = array(0, c(2, 2, 3))
#' OP[,,1] = diag(2) * exp(i)/2
#' OP[,,2] = diag(2)
#' OP[,,3] = diag(2) * exp(-i)/2
#' freq = c(-pi/3, 0, pi/3)
#' A = freqdom(OP, freq)
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