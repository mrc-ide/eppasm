# Natural age to index
a2i <- function(x, min=15, max=80) which(min:max %in% x) 

# Non number to value
na2num <- function(x, y) {x[is.na(x)] <- y; return(x)}

# new sex id for old, new, and risk groups
sid  <- function(x, fp) if (x==2) (length(fp$pi)+1):(length(fp$pi)*2) else 1:length(fp$pi)

# Make conformable arrays if there > 1 risk group
updateRiskGroup <- function(fp, mx) {
  fp$basepop <- f_sa(fp$basepop,,fp)
  if ( any(fp$cd4_mort > 1) ) fp$cd4_mort <- fp$cd4_mort/max(fp$cd4_mort)
  fp$cd4_mort <- f_sa(fp$cd4_mort,TRUE,fp)
  fp$cd4_prog <- f_sa(fp$cd4_prog,TRUE,fp)
  fp$cd4_initdist <- f_sa(fp$cd4_initdist,TRUE,fp)
  fp$art_mort <- f_sa(fp$art_mort,TRUE,fp)
  return(fp)
}

# Koronecker delta
krono <- function(x, y) ifelse(x==y, 1, 0)

# Scale to 1 by row or column
scale2one <- function(x, dim=1) sweep(x, dim, apply(x, dim, sum), '/')

# sweep with multiply function by default
sweepX <- function(...) sweep(..., FUN="*")

# Age-mixing, fraction of pns formed between indivdiual i and i'
# -----------------------------------------------------------------------------
# Log logistic 
log.logistic <- function(dis, ka = mx$kappa, rho = mx$rho, r = mx$r) 
  (ka*rho^ka)*((dis+r)^(ka-1)) / (1 + ((dis+r)^ka)^2)

# for female, male is the other way
Dkaa. <- function(a, a., mx) if (a <= a.) log.logistic(a. - a) else 0

Dmix <- function(fp, mx) {
  D <- array(0, c(fp$ss$pAG, fp$ss$pAG, 2))
  D[,,fp$ss$f.idx] <- t(sapply(fp$ss$AGE_START:fp$ss$AGE_END, 
             function(y) sapply(fp$ss$AGE_START:fp$ss$AGE_END, 
             function(x) Dkaa.(y, x, mx)))) # Female age in rows
  # Remove before first sex and after last sex
  D[,,fp$ss$f.idx][a2i(fp$ss$AGE_START:(mx$A1-1)), ] <- 0
  D[,,fp$ss$f.idx][a2i((mx$A2+1):fp$ss$AGE_END), ] <- 0
    D[,,fp$ss$f.idx] <- scale2one(D[,,fp$ss$f.idx])
    D[,,fp$ss$f.idx][is.na(D[,,fp$ss$f.idx])] <- 0  
  # Male, TODO: check for other assumptions
  D[,,fp$ss$m.idx] <- t(D[,,fp$ss$f.idx])
  return(D)
}
# mx$D <- Dmix(fp, mx)

# Risk group-mixing 
# -----------------------------------------------------------------------------
# Number of contacts for sex, age and risk group level (given as geometric mean)
Csag <- function(s, a, g, mx) with(mx,
          Mx[a2i(a), s] / tau^(sum(gamma * (seq_along(gamma)-1))) * tau^(g-1))
# Csag(f.idx, a = 15, 1, mx)

# Pop wrt sex, age, and risk group
Nsag <- function(s, a, g, pop, i, fp) sum(pop[a2i(a), sid(s, fp)[g],,i])
# Nsag(f.idx, 15, 1, pop, i)

# Mixing rate
Cmix <- function(s, a, a., g, g., pop, i, mx, fp) {
  opp <- ifelse(s == 1, 2, 1)
  num <- Csag(opp, a., g., mx) * Nsag(opp, a., g., pop, i, fp)
  den <- sum(sapply(seq_along(mx$gamma),
                    function(x) Csag(opp, a., x, mx) * Nsag(opp, a., x, pop, i, fp)))
  assort <- (1-mx$epsilon) * krono(a, a.) + mx$epsilon * num / den 
  return( Csag(s, a, g, mx) * assort * mx$D[a2i(a), a2i(a.), s] )
}
                    
# Condom use, only have ±25 and age-dependent not risk group dependent
# -----------------------------------------------------------------------------
ConAge <- function(s, a, a., mx) {
  #  more ifs to avoid big matrix
  if (a < 25 && a. < 25) out <- mx$chi[1,1]
  if (a >= 25 && a. >= 25) out <- mx$chi[2,2]
  if (a < 25 && a. >= 25) out <- ifelse(s==2, mx$chi[2,1], mx$chi[1,2])
  if (a >= 25 && a. < 25) out <- ifelse(s==2, mx$chi[1,2], mx$chi[2,1])
  return(out)
}

# Distributing base pop., hiv in, art in, base number of risk groups and %
# TODO: allow differences for female and male
# -----------------------------------------------------------------------------
f_sa <- function(inPop, mutateOnly = FALSE, fp) {
  # mutateOnly: whether to multiply with prop or just mutate the matrix
  prop <- fp$pi
  nRg <- length(prop)
  if (nRg == 1) { 
    return(inPop)
  } else {
    if (sum(prop)!=1)
      stop('Prop of risk groups do not sum to 1.')
  }
  if ( length(dim(inPop)) == 2 ) { # for basepop
    if (mutateOnly) prop <- rep(1, length(prop))
    out <- array(0, c(dim(inPop)[1], nRg*2))
    out[, 1:nRg]           <- inPop[, fp$ss$m.idx, drop=FALSE] %o% prop
    out[, (nRg+1):(nRg*2)] <- inPop[, fp$ss$f.idx, drop=FALSE] %o% prop
    return(out)
  } else if ( length(dim(inPop)) == 3 ) { # for hiv pops
    if (mutateOnly) prop <- rep(1, length(prop))
    out <- array(0, c(dim(inPop)[1:2], nRg*2))
    mal <- sapply(prop, function(x) inPop[,,fp$ss$m.idx, drop=FALSE]*x, simplify='array')
    fem <- sapply(prop, function(x) inPop[,,fp$ss$f.idx, drop=FALSE]*x, simplify='array')
    out[,,1:nRg]           <- mal
    out[,,(nRg+1):(nRg*2)] <- fem
    return(out)  
  } else if ( length(dim(inPop)) == 4 ) { # for hiv pops
    if (mutateOnly) prop <- rep(1, length(prop))
    out <- array(0, c(dim(inPop)[1:3], nRg*2))
    mal <- sapply(prop, function(x) inPop[,,,fp$ss$m.idx, drop=FALSE]*x, simplify='array')
    fem <- sapply(prop, function(x) inPop[,,,fp$ss$f.idx, drop=FALSE]*x, simplify='array')
    out[,,,1:nRg]           <- mal
    out[,,,(nRg+1):(nRg*2)] <- fem
    return(out)
  } else if ( length(inPop) == 2 ) { # for entry pop
    if (mutateOnly) prop <- rep(1, length(prop))
    out <- c(sapply(inPop, function(x) x*prop))
    return(out)
  }
}