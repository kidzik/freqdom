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
#' @example
#' Psi = 3:1 %*% t(3:1) / 20
#' X = rar(100,3,Psi)
#' matplot(X,t='l')
#'
#' DPC = dpca(X, Ndpc=1)
dpca = function(X,
                q = 30,
                Ndpc = dim(X)[2]){
  if (!is.matrix(X))
    stop("X must be a matrix")

  res = list()

  # res$mu = rep(0,dim(X)[2])
  # if (center)
  #   res$mu = colMeans(X)

  res$spec.density = spectral.density(X)
  res$filters = dpca.filters(res$spec.density, q = q, Ndpc = Ndpc)
  res$scores = dpca.scores(X, res$filters)
  res$var = dpca.var(res$spec.density)
  res$Xhat = dpca.KLexpansion(X, res$filters)
  res
}
