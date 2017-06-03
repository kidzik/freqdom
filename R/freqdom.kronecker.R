#' For given frequency-domain operators \code{F} and \code{G} (\code{\link{freqdom}}) the function \code{freqdom.kronecker} computes their Kronecker product frequency-wise.
#' 
#' Let \eqn{F = \{ F_\theta : \theta \in S \}}, \eqn{G = \{ G_\theta : \theta \in S \}},
#' where \eqn{S} is a finite grid of frequencies in \eqn{[-\pi,\pi]}, \eqn{F_\theta \in \mathbf{C}^{p_1 \times p_2}}
#' and \eqn{G_\theta \in \mathbf{C}^{q_1 \times q_2}}.
#' 
#' We define \deqn{H_\theta = F_\theta \otimes G_\theta} as a Kronecker product of \eqn{F_\theta} and \eqn{G_\theta}, i.e. \eqn{H_\theta \in \mathbf{R}^{p_1q_1\times p_2q_2}}.
#' Function \code{freqdom.kronecker} returns \eqn{H = \{ H_\theta : \theta \in S \}}.
#'  
#' @title Compute the Kronecker product of two frequency-domain operators
#' @param F a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{F_k \in \mathbf{R}^{p_1 \times p_2}}
#' on a discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @param G a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{G_k \in \mathbf{R}^{q_1 \times q_2}}
#' on a discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @return Function returns a frequency domain object (\code{\link{freqdom}}) of dimensions \eqn{L \times (p_1q_1) \times (p_2q_2)}, where \eqn{L} is the 
#' size of the evaluation grid. The elements corresponds to \eqn{F_\theta \otimes G_\theta} defined above.
#' @noRd
# @export
#' @examples
#' n = 100
#' X = rar(n,d=3)
#' Y = rar(n,d=3)
#' SX = spectral.density(X)
#' SY = spectral.density(Y)
#' R = freqdom.kronecker(SY,SX)
freqdom.kronecker = function(F,G){
  if (!is.freqdom(F))
    stop("F must be a freqdom object")
  if (!is.freqdom(G))
    stop("G must be a freqdom object")
  
  R = G
  
  D1 = dim(F$operators)
  D2 = dim(G$operators)
  D = c(0,0,0)
  D[1] = D1[1]*D2[1]
  D[2] = D1[2]*D2[2]
  D[3] = D1[3]
  R$operators = array(0,D)
  
  for (theta in 1:length(F$freq))
    R$operators[,,theta] = F$operators[,,theta] %x% G$operators[,,theta]
  R
}

oldkronprod <- `%x%`

#' #' Frequency-wise or time-wise Kronecker product. Takes two elements
#' #' \code{freqdom} or \code{timedom} and applies the Kronecker product on
#' #' each frequency or time point. If objects of other type are provided
#' #' then the standard function is applied.
#' #'  
#' #' @title Frequency-wise or component-wise Kronecker product. 
#' #' @param e1 first element
#' #' @param e2 second element
#' #' @return object of the same type as e1 but with new dimensions
#' #' @export
#' `%x%` <- function (e1,e2) {
#'   if (is.freqdom(e1) && is.freqdom(e2))
#'     freqdom.kronecker(e1,e2)
#'   else
#'     oldkronprod(e1,e2)
#' }
