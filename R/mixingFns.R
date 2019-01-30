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

# Probability of getting infected. TODO: this is assuming all stages of HIV,
# ART treatment have the same  transmission potential; so just the prev in s,
# a, m TODO: check with the code for ts time step
# -----------------------------------------------------------------------------
Psag <- function(s, a, g, pop, i, fp) {
  num <- pop[a2i(a), sid(s, fp)[g], fp$ss$hivp.idx, i]
  den <- pop[a2i(a), sid(s, fp)[g], fp$ss$hivn.idx, i] + num
  return(num/den)
}

# Condom use, only have ±25 and age-dependent not risk group dependent
# TODO: if there is data, adjust condom level use for each level of sexual act.
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

# -----------------------------------------------------------------------------
# Extracting bits from eppasm.R
# -----------------------------------------------------------------------------
# entrant populations in H+ and H-
HnIn <- function(fp, i, pregprevlag, birthslag) with(
  fp, birthslag[,i-1] * cumsurv[,i-1] * (1-entrantprev[,i] / paedsurv_lag[i-1]) + cumnetmigr[,i-1] * (1 - pregprevlag[i-1] * netmig_hivprob))
HpIn <- function(fp, i, pregprevlag, birthslag) with(
  fp, birthslag[,i-1]*fp$cumsurv[,i-1]* entrantprev[,i] + fp$cumnetmigr[,i-1]* entrantprev[,i])

# total pop x by age groups
sumByAG <- function(x, fertile=FALSE, a = fp$ss$ag.idx, b = fp$ss$p.fert.idx) {
  if (!fertile) fastmatch::ctapply(x, a, sum) 
    else fastmatch::ctapply(x, a[b], sum)
}

# sum by age groups by columns (sexes)
sumByAGs <- function(k, fertile=FALSE) apply(k, 2, sumByAG, fertile=fertile) 

# Scale mortality cd4
scale_cd4_mort_ts <- function(fp, hivpop, artpop, i) {
  if (fp$scale_cd4_mort) {
    cd4mx_scale <- hivpop[,,,i] / (hivpop[,,,i] + colSums(artpop[,,,,i]))
    cd4mx_scale[!is.finite(cd4mx_scale)] <- 1.0
    cd4_mort_ts <- fp$cd4_mort * cd4mx_scale
  } else
    cd4_mort_ts <- fp$cd4_mort
  return(cd4_mort_ts)
}

# hiv deaths at ts
calc.agdist <- function(x) {d <- x / rep(sumByAG(x), fp$ss$h.ag.span); d[is.na(d)] <- 0; d}
f_hdeaths_p.ts <- function(cd4_mort_ts, fp, pop, hivpop, artpop, i) {
  hivdeaths.ts <- fp$ss$DT * (colSums(cd4_mort_ts * hivpop[,,,i]) + colSums(fp$art_mort * fp$artmx_timerr[, i] * artpop[,,,,i],,2))
  hivdeaths_p.ts <- apply(hivdeaths.ts, 2, rep, fp$ss$h.ag.span) * 
    apply(pop[,,fp$ss$hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age
  # round(cbind(pop[,,hivp.idx,i], hivdeaths_p.ts, pop[,,hivp.idx,i] -hivdeaths_p.ts), 4)[1:10,]
  return(hivdeaths_p.ts)
}

# eligible for art
f_artcd4_percelig <- function(fp, i) with(fp, 
  1 - (1-rep(0:1, times=c(artcd4elig_idx[i]-1, ss$hDS - artcd4elig_idx[i]+1))) * (1-rep(c(0, who34percelig), c(2, ss$hDS-2))) * (1-rep(specpop_percelig[i], ss$hDS))
)

# pregnant women
f_birthdist <- function(births.by.h.age, sid, fp, i, pop, hivpop, artpop) {
  s = sid(fp$ss$f.idx, fp)
  newBornH <- sweep(hivpop[,,,i][,fp$ss$h.fert.idx, s, drop=F], 1:2, fp$frr_cd4[,,i], "*")
  nFfertileHn <- sumByAGs(pop[,,,i][fp$ss$p.fert.idx,s,fp$ss$hivn.idx,drop=F], T)
  newBornA <- colSums(sweep(artpop[,,,,i][,,fp$ss$h.fert.idx,s,drop=F], 1:3, fp$frr_art[,,,i], "*"))
  sweep(newBornH, 2:3, births.by.h.age / (nFfertileHn + colSums(newBornH) + colSums(newBornA)), '*')
}

updatePreg <- function(art15plus.elig, births.by.h.age, sid, fp, i, pop, hivpop, artpop) {
  elGrp <- 1:(fp$artcd4elig_idx[i]-1)
  elAg <- fp$ss$h.fert.idx - min(fp$ss$h.age15plus.idx) + 1 # ???
  birthdist <- f_birthdist(births.by.h.age, sid, fp, i, pop, hivpop, artpop)
  art15plus.elig[elGrp, elAg, sid(fp$ss$f.idx, fp)] %<>% +(birthdist[elGrp,,])
  return(art15plus.elig)
}

## calculate number to initiate ART based on number or percentage
## The same for risk groups
f_artInit <- function(artpop_curr_g, art15plus.elig, fp, i, ii, sid) {
  out <- numeric(fp$ss$NG)
  DT  <- fp$ss$DT
  w. <- ifelse(DT*ii < 0.5, 0, 1)
  years <- i-((2:1)-w.)
  trans <- c(1-(DT*ii+0.5-w.), DT*ii+0.5-w.)
  for(g in 1:2) {
    if(! any( fp$art15plus_isperc[g, years] ) ) {  # both number
      out[sid(g, fp)] <- sum(fp$art15plus_num[g, years] * trans)
    } else if(all(fp$art15plus_isperc[g, years])){  # both percentage
      artcov.ii <- sum(fp$art15plus_num[g, years] * trans)
      out[sid(g, fp)] <- artcov.ii * (sum(art15plus.elig[,,sid(g, fp)]) + 
        artpop_curr_g[sid(g, fp)])
    } else if(!fp$art15plus_isperc[g,i-2-w.] & fp$art15plus_isperc[g,i-1-w.]){ # transition number to percentage
      curr_coverage <- sum(artpop_curr_g[sid(g, fp)]) / (sum(art15plus.elig[,,sid(g, fp)]) + sum(artpop_curr_g[sid(g, fp)]) )
      artcov.ii <- curr_coverage + (fp$art15plus_num[g, i-1-w.] - curr_coverage) * DT/(0.5 + w. -DT*(ii-1))
      out[sid(g, fp)] <- artcov.ii * (sum(art15plus.elig[,,sid(g, fp)]) + artpop_curr_g[sid(g, fp)])
    }
  }
  return(out)
}
# f_artInit(artpop_curr_g, art15plus.elig, fp, DT, i, ii)

## calculate ART initiation distribution
f_artDist <- function(art15plus.elig, art15plus.inits, fp, i) {
  if(!fp$med_cd4init_input[i]){
    if(fp$art_alloc_method == 4L) { ## by lowest CD4
      ## Calculate proportion to be initiated in each CD4 category
      artinit <- array(0, dim(art15plus.elig))
      for(m in fp$ss$hDS:1){
        elig_hm <- colSums(art15plus.elig[m,,])
        init_prop <- ifelse(elig_hm == 0, elig_hm, pmin(1.0, art15plus.inits / elig_hm, na.rm=TRUE))
        artinit[m,,] <- sweep(art15plus.elig[m,,], 2, init_prop, "*")
      }
    } else {
      mort.w <- sweep(fp$cd4_mort, 3, colSums(art15plus.elig * fp$cd4_mort,,2), "/")
      artinit.w <- sweep(fp$art_alloc_mxweight * mort.w, 3, (1 - fp$art_alloc_mxweight) / colSums(art15plus.elig,,2), "+")
      artinit <- pmin(sweep(artinit.w * art15plus.elig, 3, art15plus.inits, "*"),
                      art15plus.elig)      
    }
  } else {
    CD4_LOW_LIM <- c(500, 350, 250, 200, 100, 50, 0)
    CD4_UPP_LIM <- c(1000, 500, 350, 250, 200, 100, 50)

    j <- fp$med_cd4init_cat[i]
    medcat_propbelow <- (fp$median_cd4init[i] - CD4_LOW_LIM[j]) / (CD4_UPP_LIM[j] - CD4_LOW_LIM[j])
    elig_below <- colSums(art15plus.elig[j,,,drop=FALSE],,2) * medcat_propbelow
    if(j < fp$ss$hDS) elig_below <- elig_below + colSums(art15plus.elig[(j+1):fp$ss$hDS,,,drop=FALSE],,2)
    elig_above <- colSums(art15plus.elig[j,,,drop=FALSE],,2) * (1.0-medcat_propbelow)
    if(j > 1) elig_above <- elig_above + colSums(art15plus.elig[1:(j-1),,,drop=FALSE],,2)

    initprob_below <- pmin(art15plus.inits * 0.5 / elig_below, 1.0, na.rm=TRUE)
    initprob_above <- pmin(art15plus.inits * 0.5 / elig_above, 1.0, na.rm=TRUE)
    initprob_medcat <- initprob_below * medcat_propbelow + initprob_above * (1-medcat_propbelow) # ??? prob > 1

    artinit <- array(0, dim=c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG))

    if(j < fp$ss$hDS)
      artinit[(j+1):fp$ss$hDS,,] <- sweep(art15plus.elig[(j+1):fp$ss$hDS,,,drop=FALSE], 3,
                                          initprob_below, "*")
    artinit[j,,] <- sweep(art15plus.elig[j,,,drop=FALSE], 3, initprob_medcat, "*")
    if(j > 0)
      artinit[1:(j-1),,] <- sweep(art15plus.elig[1:(j-1),,,drop=FALSE], 3, initprob_above, "*")
  }
  return(artinit)
}
# f_artDist(fp, art15plus.elig, art15plus.inits, i)

# End bits from eppasm.R
# -----------------------------------------------------------------------------
