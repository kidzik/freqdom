#@exportClass freqdom
#setClass("freqdom", representation(operators = "array", freq = "vector"))

#' Creates an object of class  \code{freqdom}. This object corresponds to a functional with domain \eqn{[-\pi,\pi]} and some complex vector space as codomain.
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
#' @title Create an object corresponding to a frequency domain functional
#'
#' @param F a vector, a matrix or an array. For vectors \eqn{F[k], 1\leq k\leq K} are complex numbers. For matrices \eqn{F[k,]} are complex vectors. For arrays the elements \eqn{F[,,k]}, are complex valued \eqn{(d_1\times d_2)} matrices (all of same dimension).
#' @param freq a vector of dimension \eqn{K} containing frequencies in \eqn{[-\pi,\pi]}.
#' @return Returns an object of class \code{\link{freqdom}}. An object of class  \code{\link{freqdom}} is a list containing the following components:
#' * \code{operators} \eqn{\quad} the array \code{F} as given in the argument.
#' * \code{freq} \eqn{\quad} the vector \code{freq} as given in the argument.
#' @seealso \code{\link{fourier.transform}}
#' @keywords classes
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
  res = list()
  if (is.vector(F)){
	if(length(F) != length(freq))
	stop("length of F must match number of frequencies")
    res$operators = array(0,c(1,1,length(F)))
    res$operators[1,1,] = F
  }
  else if (is.matrix(F)){
	if(length(freq) != dim(F)[1])
	stop("number of rows of F must match number of frequencies")
    res$operators = array(0,c(1,dim(F)[2],dim(F)[1]))
    res$operators[1,,] = t(F)
  }
  else if (is.array(F) && length(dim(F)) == 3)
  {
	if(length(freq) != dim(F)[3])
	stop("number of matrices of F must match number of frequencies")
    res$operators = F
  }
  else{
    stop("F must be a vector, a matrix or an array")
  }
  if(length(freq)==1)
  dimnames(res$operators)[[3]]<-list(paste("frequency", freq))
  dimnames(res$operators)[[3]]<-paste("frequency", freq)
  res$freq = freq
  class(res) = 'freqdom'
  res
}
