#' Computes the proportion and cumulative proportion of variance explained by a given number of dynamic principal components.
#' 
#' Consider a spectral density matrix \eqn{\mathcal{F}_\omega} and let \eqn{\lambda_\ell(\omega)} by the \eqn{\ell}-th dynamic eigenvalue. The proportion of variance described by the \eqn{\ell}-th dynamic principal component is given as \deqn{v_\ell:=\int_{-\pi}^\pi \lambda_\ell(\omega)d\omega/\int_{-\pi}^\pi \mathrm{tr}(\mathcal{F}_\omega)d\omega}. This function numerically computes the vectors \eqn{(v_\ell\colon 1\leq \ell\leq d)} and \eqn{\big(\sum_{j=1}^\ell v_j\colon 1\leq \ell\leq d\big)}.
#'
#' For more details we refer to Chapter 9 in Brillinger (2001), Chapter 7.8 in Shumway and Stoffer (2006) and to H\"ormann et al. (2015).
#'
#' @title Proportion of variance explained.
#' @param F \eqn{(d\times d)} spectral density matrix, provided as an object of class \code{freqdom}. To guarantee accuracy of numerical integration it is important that \eqn{F\$}freq is a dense grid of frequencies in \eqn{[-\pi,\pi]}.
#' @return A \eqn{(2\times d)}-matix. The first row contains the \eqn{v_\ell} the second row contains the cumulative sums of the \eqn{v_\ell}.
#' @seealso \code{\link{dpca.scores}}, \code{\link{dpca.filters}}, \code{\link{spectral.density}},
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R.H., and Stoffer, D.S.
#' \emph{Time Series Analysis and Its Applications} (2006), Springer, New York.
#' @export
dpca.var = function(F){
  stop("Not implemented")
}
