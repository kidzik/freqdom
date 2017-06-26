freqdom.lags = function(X){
  if (is.null(X$freq))
    lags = X$lags
  else
    lags = X$freq
  lags
}
