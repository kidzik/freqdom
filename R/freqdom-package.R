#' Frequency domain methods for stationary multivariate time series
#'
#' Methods for analyizing stationary multivariate time series, including dynamic principal components
#' and dynamic linear models. Package also provides a toolset for developers simplifying
#' construction of new frequency domain based methods for multivariate signals.
#'
#' \pkg{freqdom} package allows you to seamlessly manipulate time series objects
#' in both time and frequency domains. We implement several multivariate linear
#' methods:
#' * Dynamic principal component analysis, enabling spectral decomposition of
#' a multivariate time series into to uncorrelated components. Compared to
#' classical PCA, these components are uncorrelated for any lag. For details refer to
#' functions
#' \code{\link{dpca}} for estimating the components, \code{\link{dpca.scores}} for estimating
#' scores and \code{\link{dpca.KLexpansion}} for retrieving the signal from components,
#'
#' Apart from implementing frequency domain techniques for stationary multivariate time series,
#' \pkg{freqdom} provides a toolset of operators such as the multivariate Fourier Transform
#' (\code{\link{fourier.transform}}) or a multivariate spectral density operator
#' (\code{\link{spectral.density}}) as well as simulation of popular multivariate time series
#' models \code{\link{rar}}, \code{\link{rma}} generating multivariate
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
