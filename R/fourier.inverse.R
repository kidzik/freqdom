# get one lag: A_lag = \int R e^{itlag}
inv.fourier.one = function(R,lag){
  A = array(0,dim(R$operators)[1:2])
  A[,] = 0
  for (theta in 1:length(R$freq))
    A = A + R$operators[,,theta] * exp(R$freq[theta]*1i*lag)
  A / length(R$freq)
}

#' Computes Fourier coefficients of some functional represented by an object of class \code{freqdom}.
#'
#' Consider a function \eqn{F \colon [-\pi,\pi]\to\mathbf{C}^{d_1\times d_2}}. Its \eqn{k}-th Fourier
#' coefficient is given as
#' \deqn{
#'   \frac{1}{2\pi}\int_{-\pi}^\pi F(\omega) \exp(ik\omega)d\omega.
#' }
#' We represent the function \eqn{F} by an object of class \code{freqdom} and approximate the integral via
#' \deqn{
#' \frac{1}{|F\$freq|}\sum_{\omega\in {F\$freq}} F(\omega) \exp(i k\omega),
#' }
#' for \eqn{k\in} lags.
#'
#' @title Coefficients of a discrete Fourier transform
#' @param F an object of class \code{\link{freqdom}} which is corresponding to a function with values in \eqn{\mathbf{C}^{d_1\times d_2}}. To guarantee accuracy of inversion it is important that \code{F}\eqn{\$}\code{freq} is a dense grid of frequencies in \eqn{[-\pi,\pi]}.
#' @param lags lags of the Fourier coefficients to be computed.
#' @return An object of class \code{\link{timedom}}. The list has the following components:
#' * \code{operators} \eqn{\quad} an array. The \eqn{k}-th matrix in this array corresponds to the \eqn{k}-th Fourier coefficient.
#' * \code{lags} \eqn{\quad} the lags of the corresponding Fourier coefficients.
#' @export
#' @seealso \code{\link{fourier.transform}}, \code{\link{freqdom}}
#' @keywords time.domain frequency.domain
#' @examples
#' Y = rar(100)
#' grid = c(pi*(1:2000) / 1000 - pi) #a dense grid on -pi, pi
#' fourier.inverse(spectral.density(Y, q=2, freq=grid))
#'
#' # compare this to
#' cov.structure(Y)

fourier.inverse = function(F,lags=0){
  if (!is.freqdom(F))
    stop("F must be an object of class freqdom")
  if (!is.numeric(lags) || !is.vector(lags))
    stop("lags must be a vector of integers")

  H = length(lags)
  A = array(0, dim=c(dim(F$operators)[1:2],H))

  # TODO: this should be FFT
  for (h in 1:H)
    A[,,h] = inv.fourier.one(F,lags[h])

  timedom(Re(A),lags)
}
