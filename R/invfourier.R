# get one lag: A_lag = \int R e^{itlag}
inv.fourier.one = function(R,lag){
  A = array(0,dim(R$operators)[2:3])
  A[,] = 0
  for (theta in 1:length(R$freq))
    A = A + R$operators[theta,,] * exp(R$freq[theta]*1i*lag)
  A / length(R$freq)
}

#' For a given set of lags \eqn{G \subset \mathbf{Z}} and a set of operators
#' \eqn{R = \{ R_k : k \in S \}} (\code{\link{freqdom}}), where \eqn{S \subset [-\pi,\pi]},
#' the function \code{invfourier} evaluates the Inverse Fourier transform of
#' \eqn{R} on the set of lags \eqn{G}.
#' 
#' Given a series of lags \eqn{G \subset \mathbf{Z}} and a set of operators
#' \eqn{R = \{ R_\theta : \theta \in S \}} (\code{\link{freqdom}}), where \eqn{S \subset [-\pi,\pi]},
#' the function \code{invfourier} evaluates the Fourier transform of \eqn{R}  on each \eqn{\theta \in G}.
#' More precisely, for each \eqn{k \in G} it computes
#' \deqn{ A_k = \frac{1}{|S|} \sum_{\theta \in S} R_\theta e^{i\theta k}.}
#' 
#' It returns a time-domain operator \eqn{A = \{ A_k : k \in G \}}.
#' @title Inverse Fourier Transform of a given frequency-domain filter.
#' @param R a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{R_k \in \mathbf{R}^{p \times p}}
#' on a discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @param lags a vector lags on which the inverse Fourier transform should be evaluated. 
#' @return A time-domain object (\code{\link{timedom}}), a set of linear operators \eqn{A_k \in \mathbf{R}^{p_1 \times p_2}}
#' for selected lags \code{lags}. Every operator \eqn{A_k} is the inverse Fourier transform of \eqn{R} evaluated at lag \eqn{k}.
#' @export
#' @seealso \code{\link{fourier.transform}}
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' #estimate regressors in model $Y_t = \sum_{i\in Z} A_i X_{t-i}$
#' SYX = spectral.density(Y, X)
#' SXX = spectral.density(X)
#' R = freqdom.ratio(SYX,SXX, n)
#' A = invfourier(R) 
invfourier = function(R,lags=0:0){
  if (!is.freqdom(R))
    stop("R must be a freqdom object")
  if (!is.numeric(lags) || !is.vector(lags))
    stop("lags must be a vector of integers")
  
  H = length(lags)
  A = array(0, dim=c(H, dim(R$operators)[2:3]))

  # TODO: this should be FFT
  for (h in 1:H)
    A[h,,] = inv.fourier.one(R,lags[h])
  
  timedom(Re(A),lags)
}
