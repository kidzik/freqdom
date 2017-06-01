#' Computes the frequency response function of a linear filter.
#' 
#' Consider a filter (a sequence of vectors or matrices) \eqn{(A_k)_{k\in A\$lags}}. Then this function computes
#' \deqn{\sum_{k\in A\$lags} A_k e^{-ik\omega}}
#' for all frequencies \eqn{\omega} listed in the vector \code{freq}.
#'
#' @title Computes the Fourier transform of a filter given as \code{timedom} object.
#' @param A an object of class \code{timedom}.
#' @param freq a vector of frequencies \eqn{\in [-\pi, \pi]}. 
#' @return An object of class \code{freqdom}.
#' @seealso \code{\link{fourier.inverse}}
#' @examples
#' X = rar(100)
#' C = cov.structure(X,lags=-2:2)
#' F = fourier.transform(C) # a simple spectral density estimator
#' Cinv = fourier.inverse(F)
#' @export
fourier.transform = function(A,freq=pi*-100:100/100){
  if (is.vector(A) || is.matrix(A))
    A = timedom(A)
  if (!is.timedom(A))
    stop("A must be a time domain object")
  
  thetas = freq

  if (!is.vector(thetas) || !is.numeric(thetas) ||
        max(thetas) > pi + 0.000001 || min(thetas) < - pi - 0.000001)
    stop("freq must be a vector of real number from [-pi,pi] intercval")
  
  D = dim(A$operators)
  nbasisX = D[2]
  nbasisY = D[3]
  lags = A$lags
  
  S = array(0,c(length(thetas),nbasisX,nbasisY))
  
  # Compute sum at each frequency (TODO: this should be FFT whenever possible)
  for (theta in 1:length(thetas))
  {
    # Compute one summand
    H = length(A$lags)
    
    for (h in 1:H){
      lag = lags[h]
      R = exp(-1i*lag*thetas[theta]) * A$operators[h,,]
      S[theta,,] = S[theta,,] + R
    }
    S[theta,,] = (S[theta,,]) #/(max(thetas) - min(thetas)) # TODO: CHECK
  }
  freqdom(S,thetas)
}
