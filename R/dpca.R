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
#' @import astsa
#' @examples
#' library("astsa")
#' Y = scale(lap)
#'
#' # Compute DPCA with only one component
#' res.dpca = dpca(Y, q = 5, Ndpc = 1)
#'
#' # Compute PCA with only one component
#' res.pca = prcomp(Y,center = TRUE)
#' res.pca$x[,-1] = 0
#'
#' # Reconstruct the data
#' var.dpca = (1 - sum((res.dpca$Xhat - Y)**2) / sum(Y**2))*100
#' var.pca = (1 - sum( (res.pca$x %*% t(res.pca$rotation) - Y)**2) / sum( (Y)**2))*100
#'
#' cat("Variance explained by DPCA:\t",var.dpca,"%\n")
#' cat("Variance explained by PCA:\t",var.pca,"%\n")
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
