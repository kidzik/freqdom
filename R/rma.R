#' Generates a zero mean vector moving average process.
#' 
#' We simulate a vector moving average process
#' \deqn{
#'   X_t=\varepsilon_t+\sum_{k \in lags} \Psi_k \varepsilon_{t-k},\quad 1\leq t\leq n.
#' }
#' The innovation process \eqn{\varepsilon_t} is either multivariate normal or multivarite \eqn{t} with
#' a predefined covariance/scale matrix sigma and zero mean. The noise is generated with the
#' package \code{mvtnorm}. For Gaussian noise we use \code{\link[mvtnorm]{rmvnorm}}. For Student-t noise we use
#' \code{\link[mvtnorm]{rmvt}}. The parameters sigma and df are imported as arguments, otherwise we use default settings.
#'
#' @title Moving avarege process
#' @param n number of observations to generate.
#' @param d dimension of the time series.
#' @param Psi a \code{\link{timedom}} object with operators \code{Psi$operators}, where \code{Psi$operators[,,k]}
#' is the operator on thelag \code{lags[k]}. If no value is set then we generate a vector moving average process
#' of order \eqn{1}. Then, \code{Psi$lags = c(1)} and \code{Psi$operators[,,1]} is proportional to \eqn{\exp(-(i+j)\colon 1\leq i, j\leq d)} and such
#' that the spectral radius of \code{Psi[,,1]} is \eqn{1/2}.
#' @param noise \code{mnormal} for multivariate normal noise or \code{mt} for multivariate \eqn{t} noise. If not specified \code{mnormal} is chosen.
#' @param sigma covariance  or scale matrix of the innovations.
#' @param df degrees of freedom if \code{noise = "mt"}.
#' @return A matrix with d columns and n rows. Each row corresponds to one time point.
#' @export
#' @seealso \code{\link{rar}}
rma = function(n, d = 2, Psi = NULL, noise = c("mnormal","mt"), sigma = NULL, df = 4)
{
  if (n < 1)
    stop ("n must be a positive integer")
  if (d < 1)
    stop ("d must be a positive integer")
  if (is.null(d))
    stop("Can't determine the dimension. Specify d or give Psi.")
  if (is.null(Psi)){
  	Psi = exp(-(1:d)) %*% t(exp(-(1:d)))
  	Psi = Psi/norm(Psi,type="2")/2
  	Psi = timedom(Psi,lags=0)
  }
  if (is.null(sigma))
    sigma = diag(d)
  
  lag = rev(Psi$lags)[1] - Psi$lags[1]
  X = rar(n + lag, d=d, Psi = matrix(0,d,d), noise=noise, sigma = sigma, df=df)
  print(X)
  print(Psi)
  print(is.timedom(Psi))
  print(lag)
  linproc(X[lag + 1:n,], Psi, noise=function(d){rep(0,d)})
}

rma.old = function(n, lag=2, d=NULL, noise=NULL)
{
  rma.proc(rar(n+lag-1,d=d,noise=noise))
}

#' @noRd
# @export
rma.proc = function(TS, lag=2){
  n = dim(TS)[1]
  RES = c()
  for (x in 1:(n-lag+1)){
    new = TS[x + 0:(lag-1),]
    RES = rbind(RES,colMeans(new))
  }
  RES
}
