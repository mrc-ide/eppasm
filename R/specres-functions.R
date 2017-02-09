incid.specres <- function(x) colSums(x$newinf.m[4:10,]+x$newinf.f[4:10,]) / colSums(x$totpop.m[4:10,]+x$totpop.f[4:10,]-(x$hivnum.m[4:10,]+x$hivnum.f[4:10,]))

prev.specres <- function(x) colSums(x$hivnum.m[4:10,]+x$hivnum.f[4:10,])/colSums(x$totpop.m[4:10,]+x$totpop.f[4:10,])

aidsdeaths.specres <- function(x) colSums(x$aidsdeaths.m[-(1:3),]+x$aidsdeaths.f[-(1:3),])

agemx.specres <- function(specres){

  deaths <- with(specres, natdeaths+hivdeaths)
  pop <- with(specres, (totpop[,,-1]+totpop[,,-dim(totpop)[3]])/2)

  mx <- array(0, dim=dim(deaths), dimnames(deaths))
  mx[,,-1] <- deaths[,,-1] / pop

  return(mx)
}

natagemx.specres <- function(specres){

  deaths <- specres$natdeaths
  pop <- with(specres, (totpop[,,-1]+totpop[,,-dim(totpop)[3]])/2)

  mx <- array(0, dim=dim(deaths), dimnames(deaths))
  mx[,,-1] <- deaths[,,-1] / pop

  return(mx)
}

calc_nqx.specres <- function(specres, n=45, x=15, nonhiv=FALSE){
  if(nonhiv)
    mx <- natagemx(specres)
  else
    mx <- agemx(specres)
  1-exp(-colSums(mx[as.character(x+0:(n-1)),,]))
}


aggr_specres <- function(specreslist){
  out <- lapply(do.call(mapply, c(FUN=list, specreslist, SIMPLIFY=FALSE)), Reduce, f="+")
  class(out) <- "specres"
  return(out)
}
