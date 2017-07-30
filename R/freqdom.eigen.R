#' Gives the eigendecomposition of objects of class \code{freqdom}.
#'
#' This function makes an eigendecomposition for each of the matrices \code{F\$operator[,,k]}.
#'
#' @title Eigendecompose a frequency domain operator at each frequency
#' @param F an object of class freqdom. The matrices \code{F\$operator[,,k]} are required to be square matrices, say \eqn{d \times d}.
#' @return Returns a list. The list is containing the following components:
#' * \code{vectors} \eqn{\quad} an array containing \eqn{d} matrices. The \eqn{i}-th matrix contains in its \eqn{k}-th row the conjugate transpose eigenvector belonging to the \eqn{k}-th largest eigenvalue of \code{F\$operator[,,i]}.
#' * \code{values} \eqn{\quad} matrix containing in \eqn{k}-th column the eigenvalues of \code{F\$operator[,,k]}.
#' * \code{freq} \eqn{\quad} vector of frequencies defining the object \code{F}.
#' @importFrom graphics par plot title
#' @importFrom stats optim rnorm
# @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
# \emph{Dynamic functional principal components.} Journal of the Royal
# Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @keywords frequency.domain
#' @seealso \code{\link{freqdom}}
#' @export
freqdom.eigen = function(F){
  # TODO: It would be cleaner if this function returned two frequency domain objects
  if (!is.freqdom(F))
    stop("F must be an object of class freqdom")

  op = F$operators[,,1]
  if (dim(op)[1] != dim(op)[2])
    stop("square matrices required")

  thetas = F$freq
  E = list()
  nbasis = dim(F$operators)[2]
  T = length(thetas)

  E$freq = F$freq
  E$vectors = array(0,c(nbasis,nbasis,T))
  E$values = array(0,c(nbasis,T))
  Prev = diag(nbasis)

  for (theta in 1:T)
  {
    Eg = close.eigen(F$operators[,,theta],Prev)
    #Eg = eigen(F$operators[,,theta])
    Prev = Eg$vectors

    ## TAKE EIGEN WITHOUT(!) HEURISTIC
    #Eg = eigen(F$operators[,,theta])

    E$vectors[,,theta] = Conj(t(Eg$vectors))
    E$values[,theta] = Eg$values
  }

  E
  #  OneSideFilter(E,10,thetas)
}


# OneSideFilter = function(E,L,thetas,side="left"){
#   D = dim(E$vectors)
#   T = length(thetas)
#
#   if (side == "left") side = 1
#   else side = -1
#
#   weights = c(sqrt(1:(L)),-sqrt((L+1):1))
#   for (k in 1:1){
#
#     fr = function(param){
#       RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) #inv fourier
#       -sum(weights*abs(t(RES))^2)
# #      - sum(abs(sqrt((L+1):1)*t(RES[,L + 1:(L+1)]))^2) - sum(abs(sqrt(1:L)*t(RES[,1:L]))^2)
#
#     }
#
# start = rep(0,T)
# start = rnorm(T)
# Opt = optim(start,fr,method = "L-BFGS-B",control = list(maxit = 100, temp = 1000, trace=TRUE))
# #Opt = optim(rep(0,T),fr,method = "SANN",control = list(maxit = 30000, temp = 1000, trace=TRUE),
#     #  lower=-pi,upper=pi)
#
#     param = Opt$par
#     # plot(exp(1i*param)*1:(2*T+1),t='l')
#     ORG = t(E$vectors[,k,]) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
#     RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
#
#     E$vectors[,k,] = (exp(1i*param) * (E$vectors[,k,]))
#
# print("Weights prev i<0");
# print(sum(abs(t(ORG[,1:L]))^2))
# print("Weights prev i>=0");
# print(sum(abs(t(ORG[,L+1:(L+1)]))^2))
#
# print("Weights i<0");
# print(sum(abs(t(RES[,1:L]))^2))
# print("Weights i>=0");
# print(sum(abs(t(RES[,L+1:(L+1)]))^2))
#
# par(mfrow=c(2,2))
# plot(thetas, param, t='l', xlab="frequency", ylab="magnitude")
# title("Rotation")
# plot(0:(2*L)-L,weights,xlab="l",ylab="L_2 norm",type='h');
#   title("Rotation")
#
# plot(0:(2*L)-L,colMeans(abs(ORG[,])^2),xlab="l",ylab="L_2 norm",type='h');
# title("Original filters' magnitudes")
# plot(0:(2*L)-L,colMeans(abs(RES[,])^2),xlab="l",ylab="L_2 norm",type='h');
# title("Rotated filters' magnitudes")
# par(mfrow=c(1,1))
#
#   }
#   E
# }
