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
#' @seealso \code{\link{freqdom.inverse}}, \code{\link{freqdom.ratio}}
#' @export
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
  if (dim(F$operators)[3] != dim(G$operators)[2])
    stop("Dimensions of operators don't match")
  
  R = G
  D = c(dim(F$operators)[1],
        dim(F$operators)[2],
        dim(G$operators)[3])
  R$operators = array(0,D)
  for (theta in 1:length(F$freq))
    R$operators[theta,,] = F$operators[theta,,] %*% G$operators[theta,,]
  R
}

oldprod <- `%*%`

#' Frequency-wise or time-wise matrix product. Takes two elements
#' \code{freqdom} or \code{timedom} and multiplies them on
#' each frequency or time point. If objects of other type are provided
#' then the standard multiplication is applied.
#'  
#' @title Frequency-wise or component-wise matrix product. 
#' @param e1 first element
#' @param e2 second element
#' @return object of the same type as e1 but with new dimensions
#' @export
`%*%` <- function (e1,e2) {
  if (is.freqdom(e1) && is.freqdom(e2))
    freqdom.product(e1,e2)
  else
    oldprod(e1,e2)
}
