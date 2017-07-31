#' For given frequency-domain operators \code{F} and \code{G} (\code{\link{freqdom}}) the function \code{freqdom.kronecker} computes their matrix product frequency-wise.
#'
#' Let \eqn{F = \{ F_\theta : \theta \in S \}}, \eqn{G = \{ G_\theta : \theta \in S \}},
#' where \eqn{S} is a finite grid of frequencies in \eqn{[-\pi,\pi]}, \eqn{F_\theta \in \mathbf{C}^{p \times q}}
#' and \eqn{G_\theta \in \mathbf{C}^{q \times r}}.
#'
#' We define \deqn{H_\theta = F_\theta G_\theta} as a matrix product of \eqn{F_\theta} and \eqn{G_\theta}, i.e. \eqn{H_\theta \in \mathbf{R}^{p\times r}}.
#' Function \code{freqdom.product} returns \eqn{H = \{ H_\theta : \theta \in S \}}.
#'
#' @title Compute a matrix product of two frequency-domain operators
#' @param F frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{F_\theta \in \mathbf{R}^{p \times q}} defined
#' on a discreet grid defined \eqn{S \subset [-\pi,\pi]}.
#' @param G frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{G_\theta \in \mathbf{R}^{q \times r}} defined
#' on a discreet grid defined \eqn{S \subset [-\pi,\pi]}.
#' @return Function returns a frequency domain object (\code{\link{freqdom}}) of dimensions \eqn{L \times p \times r}, where \eqn{L} is the
#' size of the evaluation grid. The elements correspond to \eqn{F_\theta * G_\theta} defined above.
# @seealso \code{\link{freqdom.inverse}}, \code{\link{freqdom.ratio}}
#' @describeIn freqdom.product Frequency-wise matrix product of two frequency-domain operators
#' @export
#' @keywords internal
#' @examples
#' n = 100
#' X = rar(n)
#' Y = rar(n)
#' SX = spectral.density(X)
#' SY = spectral.density(Y)
#' R = freqdom.product(SY,SX)
freqdom.product = function(F,G){
  if (!is.freqdom(F))
    stop("F must be a freqdom object")
  if (!is.freqdom(G))
    stop("G must be a freqdom object")
  if (dim(F$operators)[2] != dim(G$operators)[1])
    stop("Dimensions of operators don't match")

  R = G
  D = c(dim(F$operators)[1],
        dim(G$operators)[2],
        dim(F$operators)[3])
  R$operators = array(0,D)
  for (theta in 1:length(F$freq))
    R$operators[,,theta] = F$operators[,,theta] %*% (G$operators[,,theta])
  R
}

oldprod <- `%*%`

#' Frequency-wise product of freqdom objects
#'
#' @title Frequency-wise product of freqdom objects
#' @export
#' @keywords internal
`%*%` <- function(A, B){
  if (is.freqdom(A) && is.freqdom(B))
    freqdom.product(A, B)
  else{
    oldprod(A,B)
  }
}

# @describeIn freqdom.product Convenience operator for \code{freqdom.product} function
# @title Compute a matrix product of two frequency-domain operators
# @noRd
# @export
# @exportMethod %*%
#setGeneric("%*%")
#setMethod("%*%", signature("freqdom"), function(x, y) freqdom.product(x,y))
