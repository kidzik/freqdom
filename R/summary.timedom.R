#' Print timedom object
#' 
#' @title Print timedom object
# @export print timedom
#' @export
#' @keywords internal
print.timedom = function (x, ...) {
  dots <- list(...)
  operators = TRUE
  if ("operators" %in% names(dots))
    operators = dots$operators
  
  cat("Time domain object\n\n")
  if (operators)
    for (i in 1:length(x$lags)){
    cat(paste("$operators[",i,",,] (operator at the lag ",x$lags[i],")\n",sep=""))
    print(x$operators[,,i])
    cat("\n")
  }
  cat(paste("Dimentions:",paste(dim(x$operators)[-1],collapse=" x "),"\n"))
  cat(paste("\n$lags (available lags)\n",paste(x$lags,collapse=" "),"\n"))
}

#' Print freqdom object
#' 
#' @title Print freqdom object
# @export print freqdom
#' @export
#' @keywords internal
print.freqdom = function (x, ...) {
  dots <- list(...)
  operators = TRUE
  if ("operators" %in% names(dots))
    operators = dots$operators
  
  cat("Frequency domain object\n\n")
  cat(paste("Operator dimentions:",paste(dim(x$operators)[-1],collapse=" x "),"\n"))
  cat(paste("\n$operators (matrices at each frequency)\n"))
  
  if (operators)
    for (i in 1:length(x$freq)){
    cat(paste("$operators[",i,",,] (operator at frequency ",x$freq[i],")\n",sep=""))
    print(signif(x$operators[,,i],4))
    cat("\n")
  }
  cat(paste("$freq (sample frequencies)\n"))
  print(signif(x$freq,4))
}

#' Print object summary
#' 
#' @title Print object summary
# @export summary timedom
#' @export
#' @keywords internal
summary.timedom = function (object, ...) {
  print.timedom(object, FALSE)
}

#' Print object summary
#' 
#' @title Print object summary
# @export summary freqdom
#' @export
#' @keywords internal
summary.freqdom = function (object, ...) {
  print.freqdom(object, FALSE)
}
