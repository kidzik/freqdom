#' Generates \eqn{n} observations of a \eqn{d}-dimensional moving average process with operators A, i.e.
#' \deqn{ Y_t = \sum_{i \in L_A} A_i X_t },
#' where L_A is a set of lags on which operators A_i are defined.
#'
#' @title Moving avarege process
#' @param n number of observations to generate
#' @param d number of dimensions of the process
#' @param Psi time domain object describing operators 
#' @param noise the underlying X process
#' @export
rma = function(n, d = 2, Psi = NULL, noise = c("mnormal","mt"), sigma = NULL, df = 4)
{
  if (n < 1)
    stop ("n must be a positive integer")
  if (d < 1)
    stop ("d must be a positive integer")
  if (is.null(d))
    stop("Can't determine the dimension. Specify d or give Psi.")
  if (is.null(Psi))
    Psi = timedom(diag(d),lags=0)
  if (is.null(sigma))
    sigma = diag(d)
  
  lag = rev(Psi$lags)[1] - Psi$lags[1]
  X = rar(n + lag, d=d, Psi = matrix(0,d,d), noise=noise, sigma = sigma, df=df)
  linproc(X[lag + 1:n,], Psi, noise=function(d){rep(0,d)})
}

rma.old = function(n, lag=2, d=NULL, noise=NULL)
{
  rma.proc(rar(n+lag-1,d=d,noise=noise))
}

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
