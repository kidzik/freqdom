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
#' # We compute the discrete Fourier transform (DFT) of a time series X_1,..., X_T.
#' 
#' X = rar(100)
#' DFT = fourier.transform(tdX) / sqrt(dim(X)[1])
#' @export
fourier.transform = function(A,freq=pi*-100:100/100){
  if (!is.timedom(A))
    stop("A must be an object of class timedom")
  
  thetas = freq

  if (!is.vector(thetas) || !is.numeric(thetas) ||
        max(thetas) > pi + 0.000001 || min(thetas) < - pi - 0.000001)
    stop("freq must be a vector of real numbers in [-pi,pi]")
  
  D = dim(A$operators)
  nbasisX = D[1]
  nbasisY = D[2]
  lags = A$lags
  H = length(A$lags)
  
  S = array(0,c(nbasisX, nbasisY, length(thetas)))
  
  # Compute sum at each frequency (TODO: this should be FFT whenever possible)
  for (theta in 1:length(thetas))
  {
    # Compute one summand  
    for (h in 1:H){
      R = exp(-1i*lags[h]*thetas[theta]) * A$operators[,,h]
      S[,,theta] = S[,,theta] + R
    }
    #S[,,theta] = (S[,,theta]) /(max(thetas) - min(thetas)) # TODO: CHECK
  }
  freqdom(S, thetas)
}
