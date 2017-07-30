#' Generates a zero mean vector autoregressive process of a given order.
#'
#' We simulate a vector autoregressive process
#' \deqn{
#'   X_t=\sum_{k=1}^p \Psi_k X_{t-k}+\varepsilon_t,\quad 1\leq t\leq n.
#' }
#' The innovation process \eqn{\varepsilon_t} is either multivariate normal or multivariate
#' \eqn{t} with a predefined covariance/scale matrix sigma and zero mean. The noise is generated
#' with the package \code{mvtnorm}. For Gaussian noise we use \code{\link[mvtnorm]{rmvnorm}}. For Student-t noise
#' we use \code{\link[mvtnorm]{rmvt}}. The parameters sigma and df are imported as arguments, otherwise we use default
#' settings. To initialise the process we set
#' \eqn{[X_{1-p},\ldots,X_{0}]=[\varepsilon_{1-p},\ldots,\varepsilon_{0}]}. When \code{burnin} is set
#' equal to \eqn{K} then, n\eqn{+K} observations are generated and the first \eqn{K} will be trashed.
#'
#' @title Simulate a multivariate autoregressive time series
#' @param n number of observations to generate.
#' @param d dimension of the time series.
#' @param Psi array of \eqn{p \geq 1} coefficient matrices. \code{Psi[,,k]} is the \eqn{k}-th coefficient. If no value is set then we generate a vector autoregressive process of order 1. Then, \code{Psi[,,1]} is proportional to \eqn{\exp(-(i+j)\colon 1\leq i, j\leq d)} and such that the spectral radius of \code{Psi[,,1]} is 1/2.
#' @param burnin an integer \eqn{\geq 0}. It specifies a number of initial  observations to be trashed to achieve stationarity.
#' @param noise \code{mnormal} for multivariate normal noise or \code{mt} for multivariate student t noise. If not specified \code{mnormal} is chosen.
#' @param sigma covariance  or scale matrix of the innovations. By default the identity matrix.
#' @param df degrees of freedom if \code{noise = "mt"}.
#' @importFrom graphics plot title
#' @return A matrix with \code{d} columns and \code{n} rows. Each row corresponds to one time point.
#' @seealso \code{\link{rma}}
#' @keywords simulations
#' @export
rar = function(n, d = 2, Psi = NULL, burnin = 10, noise = c('mnormal', 'mt'), sigma = NULL, df = 4)
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
		Psi = Psi/norm(Psi,type="F")/2
	}

	# build coefficients matrix (initially null)
	coef = matrix(0,n+burnin,d)

	fnoise = function() {
	  mvtnorm::rmvnorm(n = 1,mean = rep(0,d), sigma = sigma)
	}
	if (noise == 'mt'){
	  fnoise = function() {
	    mvtnorm::rmvt(n = 1, sigma = sigma, df=df)
	  }
	}

	D = dim(Psi)[3]
	if (is.na(D)){
		Psi = array(Psi, c(dim(Psi),1))
		D=1
	}

	coef[1,] = fnoise()
	# each next follow Y = Psi(X) + noise
	for (i in 2:(n+burnin)){
		last = min(D,i-1)
		for (j in 1:last){
			coef[i,] = coef[i,] + as.matrix(Psi[,,j]) %*% coef[i-j,]
		}
		coef[i,] = coef[i,] + fnoise()
	}

	# return the time series
	coef[(burnin+1):(burnin+n),]
}
