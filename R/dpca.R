#' Applies DPCA to a multivariate time series
#' 
#' Apply
#' 
#' @title Apply DPCA
#' @param X ...
#' @param lags ..
#' @param Ndpc ...
#' @param center ...
#' @param ... ...
#' @return A list
#' @export
dpca = function(X, 
                q = 3,
                thresh = NULL,
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
  res$xhat = t(t(dpca.KLexpansion(X, res$filters)) + res$mu)
  res
}
