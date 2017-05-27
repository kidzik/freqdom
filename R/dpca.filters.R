#' For a given process \code{X} eigendecompose it's spectral density
#' and use an inverse fourier transform to get coefficients of the optimal
#' filter. For details please refer to Hormann et al paper.
#'
#' @title Compute DPCA filter coefficients
#' @param X multivariate stationary time series
#' @param V correlation structure between coefficients of vectors (default diagonal)
#' @param lags requested filter coefficients
#' @param q window for spectral density estimation as in \code{\link{spectral.density}}
#' @param weights as in \code{\link{spectral.density}}
#' @param freq frequency grid to estimate on as in \code{\link{spectral.density}}
#' @return principal components series
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dpca.scores}}
#' @export
dpca.filters = function(F, L = 10){
  if (!is.freqdom(F))
    stop("F must be a freqdom operator")
  if (L < 0)
    stop("L must be a non-negative integer")

  lags = -L:L
  E = freqdom.eigen(F)
  
	d = dim(E$vectors)[2]

  XI = array(0,c(length(lags),d,d))

	for (component in 1:d)
    XI[,component,] = t(exp(-(F$freq %*% t(lags)) * 1i)) %*% E$vectors[,,component] / length(F$freq) 
  # TODO: this is just 'fourier.inverse' written in a fancy way - should be simplified

	timedom(Re(XI[length(lags):1,,]),lags)
}
