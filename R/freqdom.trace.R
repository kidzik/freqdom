#' @noRd
# @export
freqdom.trace = function(x){
  lags = freqdom.lags(x)
  res = list()
  res$values = rep(0,length(lags))
  for (i in 1:length(lags))
    res$values[i] = sum(diag(x$operators[,,i]))
  res$freq = x$freq
  res$lags = x$lags
  res
}

#' #' @export
#' tr.freqdom = freqdom.trace
#' 
#' #' @export
#' tr.timedom = freqdom.trace