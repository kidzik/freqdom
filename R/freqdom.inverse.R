#' For a given frequency-domain operator \code{S} (\code{\link{freqdom}}) the function \code{freqdom.inverse} computes the inverse of \eqn{S_\theta^{-1}} at
#' each frequency from the evaluation grid.
#' 
#' Let \eqn{S = \{ S_\theta : \theta \in G \}}, where \eqn{G} is some finite grid
#' of frequencies in \eqn{[-\pi,\pi]} and \eqn{S_\theta \in \mathbf{C}^{p \times p}}.
#' At each frequency \eqn{\theta \in G} function \code{freqdom.inverse} inverts
#' the matrix \eqn{S_\theta}, using spectral decomposition
#' \deqn{S_\theta = V_\theta \Lambda_\theta V_\theta',}
#' where \eqn{V_\theta} is orthogonal and \eqn{\Lambda_\theta} is diagonal.
#' 
#' Inverting \eqn{V_\theta \Lambda_\theta V_\theta'} requires inverting  
#' the diagonal matrix \eqn{\Lambda_\theta}, which is computationaly feasible, yet
#' in practice \eqn{\Lambda_\theta} is estimated from the sample
#' and inverting small \eqn{\lambda}s leads to high variance of the estimated \eqn{S_\theta}.
#' 
#' To mitigate this problem, we trucate \eqn{\lambda}s which are close to zero, using the heuristic
#' \eqn{\lambda_i = 0} for \eqn{i} such that \eqn{\lambda_i / \lambda_1 < 1/n}, where \eqn{n} is the
#' sample size used for the prediction. Let \eqn{\Lambda_{\theta,n}^{-1}} denote the matrix with only
#' aformentioned elements inverted and the rest set to zero.
#' 
#' Resulting matrix is 
#' \deqn{S_\theta^{-1} = V_\theta \Lambda^{-1}_\theta V_\theta'.}
#' 
#' @title Compute an inverse of a given frequency-domain operator at each frequency
#' @param S a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{S_k \in \mathbf{R}^{p \times p}}
#' on some discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @param n the number of observations used for estimation - the precision of inversion is calculated using this parameter. Either \code{n} or \code{K} must be given. 
#' @param K defines how many directions should be inverted (as in \code{\link{pseudoinverse}}). Either \code{n} or \code{K} must be given. 
#' @return The function returns a frequency-domain object (\code{\link{freqdom}}) of the same dimensions as \eqn{S}.
#' The elements of the object correspond to \eqn{S_\theta^{-1}} as defined in the Details section.
#' @seealso \code{\link{freqdom.product}}, \code{\link{freqdom.ratio}}
#' @references Hormann Siegfried, Kidzinski Lukasz and Kokoszka Piotr.
#' \emph{Estimation in functional lagged regression.} 
#' Journal of Time Series Analysis 36.4 (2015): 541-561.
# @noRd
# @export
# freqdom.inverse = function(S, n=NULL, K=NULL){
#   if (!is.freqdom(S))
#     stop("S must be a freqdom object")
#   if (!is.null(K) && !is.positiveint(K+1) && is.vector(K))
#     stop("K must be a nonnegative integer")
#   
#   # Compute number of directions to estimate
#   E = freqdom.eigen(S)
#   
#   if (is.null(K)){
#     th = 1/n
#     
#     K = sum(freqdom.inflambdas(E) / freqdom.inflambdas(E)[1] >= th)
#     debug.trace(paste("n =",n,", we have chosen K",K))
#   }
#   if (K <= 1 && !is.vector(K))
#     K = 1
# 
#   R = S
#   for (theta in 1:length(S$freq)){
#     th = 0
#     # Instead of inversing S we just take several eigendirections
#     
#     D = dim(S$operators)
#     A = matrix(S$operators[,,theta],D[1],D[2])
#     
#     curK = K
#     if (is.vector(K) && length(K) > 1)
#       curK = K[theta]
#     R$operators[,,theta] = pseudoinverse(A,curK)
#   }
#   R
# }
