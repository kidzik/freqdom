#' Collects DPCA related functions of this package. Allows to extract filters (\code{\link{dpca.filters}}), scores
#' (\code{\link{dpca.scores}}), the spectral density (\code{\link{spectral.density}}), etc. related to dynamic PCA of a
#' multivariate  time series.
#'
#' This function uses default settings for \code{spectral.density()} applied to data X.
#'
#' @title Apply DPCA to some vector time series.
#' @param X data matrix
#' @param q as in \code{\link{dpca.filters}}
#' @param Ndpc as in \code{\link{dpca.filters}}
#' @param center if TRUE then X will be centered before applying DPCA
#'
#' @return A list containing
#' * \code{scores} \eqn{\quad} DPCA scores (\code{\link{dpca.scores}})
#' * \code{filters} \eqn{\quad}  DPCA filters (\code{\link{dpca.filters}})
#' * \code{spec.density} \eqn{\quad}  spectral density of X (\code{\link{spectral.density}})
#' * \code{var} \eqn{\quad} amount of variance explained by dynamic principal components (\code{\link{dpca.var}})
#' * \code{Xhat} \eqn{\quad}  Karhunen-Loeve expansion using Ndpc dynamic principal components (\code{\link{dpca.KLexpansion}})
#'
#' @keywords DPCA
#' @export
dpca = function(X,
                q = 30,
                Ndpc = dim(X)[2],
                center = TRUE){
  if (!is.matrix(X))
    stop("X must be a matrix")

    res = list()
  res$mu = rep(0,dim(X)[2])
  if (center)
    res$mu = colMeans(X)
  res$spec.density = spectral.density(X)
  res$filters = dpca.filters(res$spec.density, q = q, Ndpc = Ndpc)
  res$scores = dpca.scores(t(t(X) - res$mu), res$filters)
  res$var = dpca.var(res$spec.density)
  res$Xhat = dpca.KLexpansion(res$scores, res$filters)
  res
}
