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

## Direct incidence input
f_infections_directincid <- function(pop, i, fp) {
  if(fp$incidpopage == 0L) # incidence for 15-49 population
    xid <- fp$ss$p.age15to49.idx
  else if(fp$incidpopage == 1L) # incidence for 15+ population
    xid <- fp$ss$p.age15plus.idx
  sexinc <- fp$incidinput[i] * c(1, fp$incrr_sex[i]) * sum(pop[xid,,fp$ss$hivn.idx,i-1]) / 
            (sum(pop[xid,fp$ss$m.idx,fp$ss$hivn.idx,i-1]) + fp$incrr_sex[i] * sum(pop[xid,fp$ss$f.idx,fp$ss$hivn.idx,i-1]))
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc / 
            (colSums(pop[xid,,fp$ss$hivn.idx,i-1] * fp$incrr_age[xid,,i]) / colSums(pop[xid,,fp$ss$hivn.idx,i-1])), "*")
  return(agesex.inc)
}

# Inside calc_infections_eppspectrum
f_art.ii <- function(h, p, DT, i, ii, pop, hivpop, artpop, fp) {
  wDT <- (1-DT*(ii-1))
  art.ii <- sum(artpop[,,h,,i])
  if(sum(hivpop[,h[1],,i]) + sum(artpop[,,h[1],,i])  > 0)
    art.ii <- art.ii - sum(pop[p[1],,fp$ss$hivp.idx,i] * colSums(artpop[,,h[1],,i],,2) / (colSums(hivpop[,h[1],,i],,1) + colSums(artpop[,,h[1],,i],,2))) * wDT
  if(sum(hivpop[,tail(h, 1)+1,,i]) + sum(artpop[,,tail(h, 1)+1,,i]) > 0)
    art.ii <- art.ii + sum(pop[tail(p,1)+1,,fp$ss$hivp.idx,i] * colSums(artpop[,,tail(h, 1)+1,,i],,2) / (colSums(hivpop[,tail(h, 1)+1,,i],,1) + colSums(artpop[,,tail(h, 1)+1,,i],,2))) * wDT
  return(art.ii)
}

f_agesex.inc <- function(inc.rate, p, i, pop, fp) {
  sexinc15to49.ts <- inc.rate * c(1, fp$incrr_sex[i]) * sum(pop[p,,fp$ss$hivn.idx,i]) / 
                     (sum(pop[p,fp$ss$m.idx,fp$ss$hivn.idx,i]) + fp$incrr_sex[i] * sum(pop[p,fp$ss$f.idx,fp$ss$hivn.idx,i]))
  sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts / (colSums(pop[p,,fp$ss$hivn.idx,i] * fp$incrr_age[p,,i])/colSums(pop[p,,fp$ss$hivn.idx,i])), "*")
}

f_hiv.ii <- function(p, hid, DT, i, ii, pop) {
  wDT     <- (1-DT*(ii-1))
  hiv.ii <- sum(pop[p,,hid,i])
  hiv.ii <- hiv.ii - sum(pop[p[1],,hid,i]) * wDT
  hiv.ii <- hiv.ii + sum(pop[tail(p,1)+1,,hid,i]) * wDT
  return(hiv.ii)
}
# End bits from eppasm.R
# -----------------------------------------------------------------------------