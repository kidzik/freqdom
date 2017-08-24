#' Dynamic principal component analysis (DPCA) decomposes multivariate time series into uncorrelated components. Compared
#' to classical principal components, DPCA decomposition outputs components which are uncorrelated in time, allowing
#' simpler modeling of the processes and maximizing long run variance of the projection.
#'
#' This convenience function applies the DPCA methodology and returns filters (\code{\link{dpca.filters}}), scores
#' (\code{\link{dpca.scores}}), the spectral density (\code{\link{spectral.density}}), variances (\code{\link{dpca.var}}) and
#' Karhunen-Leove expansion (\code{\link{dpca.KLexpansion}}).
#'
#' See the example for understanding usage, and help pages for details on individual functions.
#'
#' @title Compute Dynamic Principal Components and dynamic Karhunen Loeve extepansion
#'
#' @param X a vector time series given as a \eqn{(T\times d)}-matix. Each row corresponds to a timepoint.
#' @param q window size for the kernel estimator, i.e. a positive integer.
#' @param freq a vector containing frequencies in \eqn{[-\pi, \pi]} on which the spectral density should be evaluated.
#' @param Ndpc is the number of principal component filters to compute as in \code{\link{dpca.filters}}
#'
#' @references Hormann, S., Kidzinski, L., and Hallin, M.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Brillinger, D.
#' \emph{Time Series} (2001), SIAM, San Francisco.
#' @references Shumway, R., and Stoffer, D.
#' \emph{Time series analysis and its applications: with R examples} (2010), Springer Science & Business Media
#'
#' @return A list containing
#' * \code{scores} \eqn{\quad} DPCA scores (\code{\link{dpca.scores}})
#' * \code{filters} \eqn{\quad}  DPCA filters (\code{\link{dpca.filters}})
#' * \code{spec.density} \eqn{\quad}  spectral density of \code{X} (\code{\link{spectral.density}})
#' * \code{var} \eqn{\quad} amount of variance explained by dynamic principal components (\code{\link{dpca.var}})
#' * \code{Xhat} \eqn{\quad}  Karhunen-Loeve expansion using \code{Ndpc} dynamic principal components (\code{\link{dpca.KLexpansion}})
#'
#' @keywords DPCA
#' @export
#' @examples
#' X = rar(100,3)
#'
#' # Compute DPCA with only one component
#' res.dpca = dpca(X, q = 5, Ndpc = 1)
#'
#' # Compute PCA with only one component
#' res.pca = prcomp(X, center = TRUE)
#' res.pca$x[,-1] = 0
#'
#' # Reconstruct the data
#' var.dpca = (1 - sum( (res.dpca$Xhat - X)**2 ) / sum(X**2))*100
#' var.pca = (1 - sum( (res.pca$x %*% t(res.pca$rotation) - X)**2 ) / sum(X**2))*100
#'
#' cat("Variance explained by DPCA:\t",var.dpca,"%\n")
#' cat("Variance explained by PCA:\t",var.pca,"%\n")
dpca = function(X,
                q = 30,
                freq = (-1000:1000/1000)*pi,
                Ndpc = dim(X)[2]){
  if (!is.matrix(X))
    stop("X must be a matrix")

  res = list()
  res$spec.density = spectral.density(X, freq = freq)
  res$filters = dpca.filters(res$spec.density, q = q, Ndpc = Ndpc)
  res$scores = dpca.scores(X, res$filters)
  res$var = dpca.var(res$spec.density)
  res$Xhat = dpca.KLexpansion(X, res$filters)
  res
}
