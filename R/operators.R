#' Frequency-wise sum of freqdom objects
#' 
#' @title Frequency-wise sum of freqdom objects
#' @export
#' @keywords internal
plus.freqdom = function(e1,e2){
  R = e1
  lags = freqdom.lags(e1)
  for (i in 1:length(lags))
    R$operators[,,i] = e1$operators[,,i] + e2$operators[,,i]
  R
}

plus.timedom = function(e1,e2){
  R = e1
  lags = freqdom.lags(e1)
  for (i in 1:length(lags))
    R$operators[,,i] = e1$operators[,,i] + e2$operators[,,i]
  R
}

#' Frequency-wise sum of freqdom objects
#' 
#' @title Frequency-wise sum of freqdom objects
#' @export
#' @keywords internal
"+.freqdom" = function (e1,e2) plus.freqdom(e1,e2)

#' Time-wise sum of freqdom objects
#' 
#' @title Time-wise sum of freqdom objects
#' @export
#' @keywords internal
"+.timedom" = function (e1,e2) plus.freqdom(e1,e2)

minus.freqdom = function(e1,e2){
  R = e1
  R$lags = union(e1$lags,e2$lags)
  R$operators = array(0,c(dim(R$operators)[1:2],length(R$lags)))
  
  for (i in 1:length(R$lags)){
    lag = R$lags[i]
    i1 = which(e1$lags == lag)
    i2 = which(e2$lags == lag)
    if (sum(i1)==0)
      R$operators[,,i] = - e2$operators[,,i2]
    else if (sum(i2)==0)
      R$operators[,,i] = e1$operators[,,i1]
    else
      R$operators[,,i] = e1$operators[,,i1] - e2$operators[,,i2]
  }
  R
}

#' Frequency-wise difference of freqdom objects
#' 
#' @title Frequency-wise difference of freqdom objects
#' @export
#' @keywords internal
"-.freqdom" = function (e1,e2) minus.freqdom(e1,e2)

#' Time-wise sum of freqdom objects
#' 
#' @title Time-wise difference of freqdom objects
#' @export
#' @keywords internal
"-.timedom" = function (e1,e2) minus.freqdom(e1,e2)
