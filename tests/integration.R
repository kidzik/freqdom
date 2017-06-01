OP = array(0,c(2,2,3))
OP[,,1] = diag(2)
OP[,,2] = diag(2)
OP[,,3] = diag(2)

A = timedom(OP,c(-1,0,2))
A = timedom.trunc(A, 0)

B = timedom(array(diag(2),c(2,2,1)),0)
if (!all(B$operators == A$operators))
  stop("truncation error")

OP = array(0,c(2,1,3))
OP[,,1] = 1:2
OP[,,2] = 2:3
OP[,,3] = c(-5,5)

if (sum((fourier.inverse(fourier.transform(timedom(OP,-1:1)),lags = -1:1)$operators - OP)^2)> 0.01)
  stop("Fourier transform error")
