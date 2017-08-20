#' Frequency domain basde analysis: dynamic PCA
#'
#' Implementation of dynamic principle component analysis (DPCA),
#' simulation of VAR and VMA processes and frequency domain tools.
#' The package also provides a toolset for developers simplifying
#' construction of new frequency domain based methods for multivariate signals.
#'
#' \pkg{freqdom} package allows you to manipulate time series objects
#' in both time and frequency domains. We implement dynamic principal component analysis methods,
#' enabling spectral decomposition of a stationary vector time series into uncorrelated components.
#'
#' Dynamic principal component analysis enables estimation of temporal filters which
#' transform a vector time series into another vector time series with uncorrelated components,
#' maximizing the long run variance explained.
#' There are two key differnces between classical PCA and dynamic PCA:
#' * Components returned by the dynamic procedure are uncorrelated in time, i.e. for any \eqn{i \neq j}
#' and \eqn{l \in Z}, \eqn{Y_i(t)} and \eqn{Y_j(t_l)} are uncorrelated,
#' * The mapping maximizes the long run variance, which, in case of stationary vector time series, means
#' that the process reconstructed from and \eqn{d > 0} first dynamic principal components
#' better approximates your vector time series process than the first \eqn{d} classic principal components.
#'
#' For details, please refer to literature below and to help pages of functions \code{\link{dpca}}
#' for estimating the components, \code{\link{dpca.scores}} for estimating scores and
#' \code{\link{dpca.KLexpansion}} for retrieving the signal from components.
#'
#' Apart from frequency domain techniques for stationary vector time series,
#' \pkg{freqdom} provides a toolset of operators such as the vector Fourier Transform
#' (\code{\link{fourier.transform}}) or a vector spectral density operator
#' (\code{\link{spectral.density}}) as well as simulation of vector time series
#' models \code{\link{rar}}, \code{\link{rma}} generating vector
#' autoregressive and moving average respectively.
#' These functions enable developing new techniques based on the Frequency domain analysis.
#'
#' @references Hormann Siegfried, Kidzinski Lukasz and Hallin Marc.
#' \emph{Dynamic functional principal components.} Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology) 77.2 (2015): 319-348.
#' @references Hormann Siegfried, Kidzinski Lukasz and Kokoszka Piotr.
#' \emph{Estimation in functional lagged regression.}
#' Journal of Time Series Analysis 36.4 (2015): 541-561.
#' @references Hormann Siegfried and Kidzinski Lukasz.
#' \emph{A note on estimation in Hilbertian linear models.}
#' Scandinavian journal of statistics 42.1 (2015): 43-62.
#' @aliases freqdom-package
"_PACKAGE"
