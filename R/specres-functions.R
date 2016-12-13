incid.specres <- function(x) colSums(x$newinf.m[4:10,]+x$newinf.f[4:10,]) / colSums(x$totpop.m[4:10,]+x$totpop.f[4:10,]-(x$hivnum.m[4:10,]+x$hivnum.f[4:10,]))

prev.specres <- function(x) colSums(x$hivnum.m[4:10,]+x$hivnum.f[4:10,])/colSums(x$totpop.m[4:10,]+x$totpop.f[4:10,])

aidsdeaths.specres <- function(x) colSums(x$aidsdeaths.m[-(1:3),]+x$aidsdeaths.f[-(1:3),])
