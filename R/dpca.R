#' Applies DPCA to a multivariate time series
#' 
#' Apply DPCA
#' 
#' @title Apply DPCA
#' @param X data matrix
#' @param q as in \code{\link{dpca.filters}}
#' @param thresh as in \code{\link{dpca.filters}}
#' @param Ndpc as in \code{\link{dpca.filters}}
#' @param center should X be centered before applying DPCA
#' @param ...  as in \code{\link{spectral.density}}
#' 
#' @return A list containing
#' * \code{scores} DPCA scores (\code{\link{dpca.scores}}) 
#' * \code{filters} DPCA filters (\code{\link{dpca.filters}}) 
#' * \code{spec.density} spectral density of X (\code{\link{spectral.density}}) 
#' * \code{var} long-run variance explained by components of filters (\code{\link{dpca.var}}) 
#' * \code{Xhat} Karhunen-Loeve expansion using Ndpc dynamic principal components (\code{\link{dpca.KLexpansion}})
#' 
#' @export
dpca = function(X, 
                q = NULL,
                thresh = 0.9,
                Ndpc = dim(X)[2],
                center = TRUE,
                ...){
  if (length(lags) < 0)
    stop("lags must be a vector of positive length")
  if (!is.matrix(X))
    stop("X must be a matrix")

    res = list()
  res$mu = rep(0,dim(X)[2])
  if (center)
    res$mu = colMeans(X)
  res$spec.density = spectral.density(X, ...)
  res$filters = dpca.filters(res$spec.density, q = q, thresh = thresh, Ndpc = Ndpc)
  res$scores = dpca.scores(t(t(X) - res$mu), res$filters)
  res$var = dpca.var(res$spec.density)
  res$Xhat = t(t(dpca.KLexpansion(X, res$filters)) + res$mu)
  res
}
