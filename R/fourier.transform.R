#' For a given set of frequencies \eqn{G \subset [-\pi,\pi]} and a  series of operators
#' \eqn{A = \{ A_k : k \in S \}} (\code{\link{timedom}}), where \eqn{S} is a set of lags,
#' the function \code{fourier.transform} evaluates the Fourier transform of
#' \eqn{A} on the set \eqn{G}.
#' 
#' Given a series of frequencies \eqn{G \subset [-\pi,\pi]} and a series of operators
#' \deqn{(A(h): h \in \{ -q,...,0,...,q\}),}
#' where \eqn{A(h) \in \mathbf{R}^{p_1 \times p_2}},
#' the function \code{fourier.transform} evaluates the Fourier transform on each \eqn{\theta \in G}.
#' More precisely, for each \eqn{\theta \in G} it computes
#' \deqn{ F_\theta = \sum_{h=-q}^q A(h) e^{-i\theta h}.}
#' 
#' It returns a frequency-domain operator \eqn{F = \{ F_\theta : \theta \in G \}}.
#'
#' @title Computes the Fourier transform of a given series of operators or a given multivariate time series
#' @param A a time-domain object \code{\link{timedom}}, i.e. a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
#' for some set \eqn{G \subset \mathbf{Z}} of lags. Alternatively, a multivariate time series.
#' @param freq frequencies on which the transfom should be evaluated. A vector of increasing values in \eqn{[-\pi,\pi]}. 
#' If the parameter is \code{NULL}, the default grid of \eqn{201} equidistributed points on \eqn{[-\pi,\pi]} is taken
#' as the grid. Here we refer to this set as \eqn{G}.
#' @return A frequency domain operator \eqn{F = \{ F_\theta : \theta \in G \}} (\code{\link{freqdom}}) defined on the given set of frequencies \eqn{G}.
#' Every \eqn{F_\theta} corresponds to the Fourier transform of \eqn{A} evaluated on \eqn{\theta \in G}.
#' @seealso \code{\link{invfourier}}
#' @examples
#' X = rar(100,d=2)
#' C = cov.structure(X)
#' F = fourier.transform(C) # a simple spectral density estimator
#' Cinv = invfourier(F)
#' @export
fourier.transform = function(A,freq=NULL){
  if (is.vector(A) || is.matrix(A))
    A = timedom(A)
  if (!is.timedom(A))
    stop("A must be a time domain object")
  
  if (is.null(freq))
    freq = pi*-100:100/100
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
