library("freqdom")
Psi = 3:1 %*% t(3:1) / 20
X = rar(100,3,Psi)
Ndpc = 1
q = 30

res = list()
res$spec.density = spectral.density(X)
res$filters = dpca.filters(res$spec.density, q = q, Ndpc = Ndpc)
res$scores = dpca.scores(X, res$filters)
res$var = dpca.var(res$spec.density)
res$Xhat = dpca.KLexpansion(X, res$filters)
