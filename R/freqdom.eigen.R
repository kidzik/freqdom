#' For a given frequency domain operator, compute its eigenvalues and eigendirections
#' such that for close frequencies eigendirection matrices are close to each other
#' ('quasi-continuity').
#' 
#' Let \eqn{S = \{ S_\theta : \theta \in G \}}, where \eqn{G} is some finite grid
#' of frequencies in \eqn{[-\pi,\pi]} and \eqn{S_\theta \in \mathbf{C}^{p_1 \times p_2}}.
#' At each frequency \eqn{\theta \in G} function \code{freqdom.eigen} eigendecomposes
#' the matrix \eqn{S_\theta}
#' \deqn{S_\theta = V_\theta \Lambda_\theta V_\theta',}
#' where \eqn{V_\theta} is orthogonal and \eqn{\Lambda_\theta} is diagonal.
#' 
#' Note that if \eqn{v} is an eigenvector of \eqn{S_\theta}, so is \eqn{-v}.
#' Thus, the solution is not unique and in particular for 
#' \eqn{A = diag(a_1,a_2,...,_p), a_i \in\{-1,1\}} and \eqn{W_{\theta} = S_\theta A}, we have
#' \deqn{S_\theta = W_{\theta} \Lambda_\theta W_{\theta}'}
#' 
#' One can show that this ambiguity can lead to noise solutions when inverse Fourier
#' transform is applied \eqn{A}. To mitigate this problem we follow the algorithm:
#' Let \eqn{G = \{\theta_0, \theta_1, ..., \theta_d\}} and \eqn{\theta_i < \theta_j} for \eqn{i < j}.
#' We allow arbitrary \eqn{S_{\theta_0}}. Next, for every \eqn{i > 0} we choose \eqn{A}
#' such that \eqn{S_{\theta_{i-1}}} and \eqn{S_{\theta_{i}}} are close to each other or,
#' in other words, such that
#' \deqn{\|I - S_{\theta_{i-1}} S_{\theta_i} A\|_F,}
#' where \eqn{\| \cdot \|_F} denotes the Frobenius norm.
#'
#' @title Eigendevompose a frequency domain operator at each frequency
#' @param S a frequency-domain filter of type \code{\link{freqdom}}, i.e. a set of linear operators \eqn{S_k \in \mathbf{R}^{p \times p}}
#' on a discreet grid defined of \eqn{[-\pi,\pi]}. 
#' @return A a list with elements
#' * \code{$vectors} - an array of dimensions \eqn{L \times p \times p}, where \eqn{L} is the size of the grid and \code{$vectors[i,,]}
#' correpsonds to the matrix of \eqn{\mathbf{C}^{p \times p}} eigenvectors of \eqn{S_{\theta_i}}
#' * \code{$values} - a matrix of dimensions \eqn{L \times p} where \code{$values[i,]} is a vector \eqn{\mathbf{C}^{p}} of \eqn{p} eigenvalues.
#' @importFrom graphics par plot title
#' @importFrom stats optim rnorm
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @export
freqdom.eigen = function(S){
  # TODO: It would be cleaner if this function returned two frequency domain objects
  if (!is.freqdom(S))
    stop("S must be a freqdom object")
  
  op = S$operators[1,,]
  if (dim(op)[1] != dim(op)[2])
    stop("S$operators[,,theta] must be a square matrix")
  
  thetas = S$freq
  E = list()
  nbasis = dim(S$operators)[2]
  T = length(thetas)
  
  E$freq = S$freq
  E$vectors = array(0,c(T,nbasis,nbasis))
  E$values = array(0,c(T,nbasis))
  Prev = diag(nbasis)
  
  for (theta in 1:T)
  {
    Eg = close.eigen(S$operators[theta,,],Prev)
    Prev = Eg$vectors
    
    ## TAKE EIGEN WITHOUT(!) HEURISTIC
    #Eg = eigen(S$operators[theta,,])
    
    E$vectors[theta,,] = Eg$vectors
    E$values[theta,] = Eg$values
  }
  
  E
#  OneSideFilter(E,10,thetas)
}


OneSideFilter = function(E,L,thetas,side="left"){
  D = dim(E$vectors)
  T = length(thetas)
  
  if (side == "left") side = 1
  else side = -1
  
  weights = c(sqrt(1:(L)),-sqrt((L+1):1))
  for (k in 1:1){
    
    fr = function(param){
      RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) #inv fourier
      -sum(weights*abs(t(RES))^2)
#      - sum(abs(sqrt((L+1):1)*t(RES[,L + 1:(L+1)]))^2) - sum(abs(sqrt(1:L)*t(RES[,1:L]))^2)
      
    }

start = rep(0,T)
start = rnorm(T)
Opt = optim(start,fr,method = "L-BFGS-B",control = list(maxit = 100, temp = 1000, trace=TRUE))
#Opt = optim(rep(0,T),fr,method = "SANN",control = list(maxit = 30000, temp = 1000, trace=TRUE),
    #  lower=-pi,upper=pi)
    
    param = Opt$par
    # plot(exp(1i*param)*1:(2*T+1),t='l')
    ORG = t(E$vectors[,k,]) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
    RES = t(exp(1i*param) * (E$vectors[,k,])) %*% exp(-(thetas %*% t(0:(2*L)-L)) * 1i) / T
    
    E$vectors[,k,] = (exp(1i*param) * (E$vectors[,k,]))
    
print("Weights prev i<0");
print(sum(abs(t(ORG[,1:L]))^2))
print("Weights prev i>=0");
print(sum(abs(t(ORG[,L+1:(L+1)]))^2))
    
print("Weights i<0");
print(sum(abs(t(RES[,1:L]))^2))
print("Weights i>=0");
print(sum(abs(t(RES[,L+1:(L+1)]))^2))
    
par(mfrow=c(2,2))
plot(thetas, param, t='l', xlab="frequency", ylab="magnitude")
title("Rotation")
plot(0:(2*L)-L,weights,xlab="l",ylab="L_2 norm",type='h');
  title("Rotation")
    
plot(0:(2*L)-L,colMeans(abs(ORG[,])^2),xlab="l",ylab="L_2 norm",type='h');
title("Original filters' magnitudes")
plot(0:(2*L)-L,colMeans(abs(RES[,])^2),xlab="l",ylab="L_2 norm",type='h');
title("Rotated filters' magnitudes")
par(mfrow=c(1,1))

  }
  E
}
