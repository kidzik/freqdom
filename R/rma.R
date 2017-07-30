#' Generates a zero mean vector moving average process.
#'
#' This simulates a vector moving average process
#' \deqn{
#'   X_t=\varepsilon_t+\sum_{k \in lags} \Psi_k \varepsilon_{t-k},\quad 1\leq t\leq n.
#' }
#' The innovation process \eqn{\varepsilon_t} is either multivariate normal or multivarite \eqn{t} with
#' a predefined covariance/scale matrix sigma and zero mean. The noise is generated with the
#' package \code{mvtnorm}. For Gaussian noise we use \code{\link[mvtnorm]{rmvnorm}}. For Student-t noise we use
#' \code{\link[mvtnorm]{rmvt}}. The parameters \code{sigma} and \code{df} are imported as arguments, otherwise we use default settings.
#'
#' @title Moving average process
#' @param n number of observations to generate.
#' @param d dimension of the time series.
#' @param Psi a \code{\link{timedom}} object with operators \code{Psi$operators}, where \code{Psi$operators[,,k]}
#' is the operator on thelag \code{lags[k]}. If no value is set then we generate a vector moving average process
#' of order \eqn{1}. Then, \code{Psi$lags = c(1)} and \code{Psi$operators[,,1]} is proportional to \eqn{\exp(-(i+j)\colon 1\leq i, j\leq d)} and such
#' that the spectral radius of \code{Psi[,,1]} is \eqn{1/2}.
#' @param noise \code{mnormal} for multivariate normal noise or \code{mt} for multivariate \eqn{t} noise. If not specified \code{mnormal} is chosen.
#' @param sigma covariance  or scale matrix of the innovations. If NULL then the identity matrix is used.
#' @param df degrees of freedom if \code{noise = "mt"}.
#' @return A matrix with \code{d} columns and \code{n} rows. Each row corresponds to one time point.
#' @export
#' @keywords simulations
#' @seealso \code{\link{rar}}
rma = function(n, d = 2, Psi = NULL, noise = c("mnormal","mt"), sigma = NULL, df = 4)
{
	if (!is.null(Psi))
		d = dim(Psi)[1]

	if (n < 1)
	  stop ("n must be a positive integer")
	if (d < 1)
	  stop ("d must be a positive integer")
	if (length(noise)>1)
	  noise = 'mnormal'
	if (is.null(sigma))
	  sigma = diag(d)

	if (det(sigma)<0 || !isSymmetric(sigma))
		stop("sigma is not a covariance matrix.")

	# if no operator then make some default
	if (is.null(Psi)){
		Psi = exp(-(1:d))%*%t(exp(-(1:d)))
		Psi = Psi/norm(Psi,type="2")/2
	}

	if(is.array(Psi)) nlags = dim(Psi)[3]
	if(is.matrix(Psi)) nlags = 1
	Psi.new=array(0,c(dim(Psi)[1],dim(Psi)[1],nlags+1))
	Psi.new[,,1]=diag(d)
	if(nlags==1) {
	Psi.new[,,2]=Psi
	} else
	for(i in 1:nlags){ Psi.new[,,i+1]=Psi[,,i]}

	Psi.new=timedom(Psi.new,lags=0:nlags)

	innovs =  mvtnorm::rmvnorm(n = n+nlags,mean = rep(0,d), sigma = sigma)
	if (noise == 'mt'){
	innovs = mvtnorm::rmvt(n = n+nlags, sigma = sigma, df=df)
	  }

  out=filter.process(innovs,Psi.new)
  out=out[(nlags+1):(n+nlags),]
  out
}
