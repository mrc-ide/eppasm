# -----------------------------------------------------------------------------
# Extracting bits from eppasm.R
# Most are just a wrapper
# -----------------------------------------------------------------------------
# entrant populations in H+ and H-
HnIn0 <- function(fp, i, birthslag) with(
  fp, birthslag[,i-1] * cumsurv[,i-1] / paedsurv_lag[i-1] + cumnetmigr[,i-1])
HnIn <- function(fp, i, pregprevlag, birthslag) {
  with(fp, birthslag[,i-1] * cumsurv[,i-1] * (1-entrantprev[,i] / paedsurv_lag[i-1]) + cumnetmigr[,i-1] * (1 - pregprevlag[i-1] * netmig_hivprob))
}
HpIn <- function(fp, i, pregprevlag, birthslag) {
  with(fp, birthslag[,i-1]*fp$cumsurv[,i-1]* entrantprev[,i] + fp$cumnetmigr[,i-1]* entrantprev[,i])
}

#' pop x by age groups
#' 
#' When x is a vector
#' @param x population matrix/data frame
#' @param fp fix_par
#' @param fertile restrict to fertile age-group
sumByAG <- function(x, fp, fertile=FALSE) {
  if (!fertile) fastmatch::ctapply(x, fp$ss$ag.idx, sum) 
    else fastmatch::ctapply(x, fp$ss$ag.idx[fp$ss$p.fert.idx], sum)
}

#' pop x by age groups
#' 
#' When x is a data.frame/matrix
#' @param k population matrix/data frame
#' @param fix_par fix_par
#' @param fertile restrict to fertile age-group
sumByAGs <- function(k, fix_par, fertile=FALSE) apply(k, 2, sumByAG, fp=fix_par, fertile=fertile) 

# Scale mortality cd4
scale_cd4_mort <- function(fp, hivpop, artpop, i) {
  if (fp$scale_cd4_mort) {
    cd4mx <- hivpop$get(YEAR=i) / (hivpop$get(YEAR=i) + colSums(artpop$get(YEAR=i)))
    cd4mx[!is.finite(cd4mx)] <- 1.0
    cd4_mort_ts <- cd4mx * fp$cd4_mort
  } else
    cd4_mort_ts <- fp$cd4_mort
  return(cd4_mort_ts)
}

# hiv deaths at ts
calc.agdist <- function(x, fp) {
  d <- x / rep(sumByAG(x, fp), fp$ss$h.ag.span) # percent of each age/age group
  d[is.na(d)] <- 0 
  d
}

f_hiv_death <- function(MODEL, cd4_mx, cd4_mx_db, pop, i, fp,
                        hiv_pop, art_pop, hiv_db, art_db) {
  # death by age group
  artD <- art_pop$get(YEAR = i) * fp$art_mort * fp$artmx_timerr[, i]
  hivD <- cd4_mx * hiv_pop$get(YEAR = i)
  if (MODEL==2) { # add inactive deaths
    artD <- artD + art_db$get(YEAR = i) * fp$art_mort * fp$artmx_timerr[, i]
    hivD <- hivD + cd4_mx_db * hiv_db$get(YEAR = i)
  }
  dbyAG <- fp$ss$DT * (colSums(hivD) + colSums(artD,,2))
  # deaths by single-year
  pA   <- apply(pop$data[,,fp$ss$hivp.idx,i], 2, calc.agdist, fp=fp)
  dbyA <- apply(dbyAG, 2, rep, fp$ss$h.ag.span) * pA
  return(dbyA)
}

# eligible for art
f_artcd4_percelig <- function(fp, i) with(fp, 
  1 - (1 - rep(0:1, times=c(artcd4elig_idx[i]-1, ss$hDS - artcd4elig_idx[i]+1))) * (1 - rep(c(0, who34percelig), c(2, ss$hDS-2))) * (1-rep(specpop_percelig[i], ss$hDS))
)

updatePreg <- function(art_elig, birth_agrp, fp, year, pop, hivpop, artpop) {
  list2env(fp$ss, environment())
  elDS <- 1:(fp$artcd4elig_idx[year] - 1)
  elAG <- h.fert.idx - min(h.age15plus.idx) + 1

  hivp <- hivpop$data[, h.fert.idx, f.idx, year] * fp$frr_cd4[,, year]
  hivn <- sumByAG(pop[p.fert.idx, f.idx, hivn.idx, year], fp, TRUE)
  art  <- colSums(artpop$data[,,h.fert.idx, f.idx, year] * fp$frr_art[,,,year],,2)
  birthdist <- sweep(hivp, 2, birth_agrp / (hivn + colSums(hivp) + art), "*")
  art_elig[elDS, elAG, f.idx] <- art_elig[elDS, elAG, f.idx] + birthdist[elDS,]
  return(art_elig)
}

## calculate number to initiate ART based on number or percentage INPUTS
f_artInit <- function(art_curr, art_elig, fp, year, ts) {
  out    <- c(0,0)
  DT     <- fp$ss$DT
  year_w <- ifelse(DT * ts < 0.5, 0, 1)
  trans  <- DT * ts + 0.5 - year_w
  years  <- year - ((2:1) - year_w)
  for(sex in 1:2) {
    if( !any(fp$art15plus_isperc[sex, years]) ) {
      # both number
      out[sex] <- sum(fp$art15plus_num[sex, years] %*% c(1 - trans, trans))
    } else if (all(fp$art15plus_isperc[sex, years])) {
      # both percentage
      artcov <- sum(fp$art15plus_num[sex, years] %*% c(1 - trans, trans))
      out[sex] <- artcov * (sum(art_elig[,,sex]) + art_curr[sex])
    } else if (!fp$art15plus_isperc[sex, year - (2 - year_w)] & 
                fp$art15plus_isperc[sex, year - (1 - year_w)]) {
      # transition number to percentage
      curr_cov <- art_curr[sex] / (sum(art_elig[,,sex]) + art_curr[sex])
      diff_cov <- fp$art15plus_num[sex, year - (1-year_w)] - curr_cov
      artcov <- curr_cov + diff_cov * DT / (0.5 + year_w -DT * (ts - 1))
      out[sex] <- artcov * (sum(art_elig[,,sex]) + art_curr[sex])
    }
  }
  return(out)
}
# f_artInit(artpop_curr_g, art15plus.elig, fp, DT, i, ii)

## calculate ART initiation distribution
f_artDist <- function(art_elig, art_need, fp, year) {
  list2env(fp$ss, environment())
  if (!fp$med_cd4init_input[year]) {
    if (fp$art_alloc_method == 4L) { ## by lowest CD4
      ## Calculate proportion to be initiated in each CD4 category
      out <- array(0, dim(art_elig))
      for(m in fp$ss$hDS:1){
        elig_hm <- colSums(art_elig[m,,])
        init_pr <- if (elig_hm == 0) elig_hm 
                      else pmin(1.0, art_need / elig_hm, na.rm = TRUE)
        out[m,,] <- sweep(art_elig[m,,], 2, init_pr, "*")
      }
    } else { # Spectrum Manual p168--p169, 
      expect.mort.w <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
                    colSums(art_elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")
      init.w <- sweep(fp$art_alloc_mxweight * expect.mort.w, 3,
                      (1 - fp$art_alloc_mxweight) / colSums(art_elig,,2), "+")
      out <- pmin(sweep(init.w * art_elig, 3, art_need, "*"), art_elig)
    }
  } else {

    CD4_LO <- c(500, 350, 250, 200, 100, 50, 0)
    CD4_UP <- c(1000, 500, 350, 250, 200, 100, 50)
    
    j <- fp$med_cd4init_cat[year]

    pr_below <- (fp$median_cd4init[year] - CD4_LO[j]) / (CD4_UP[j] - CD4_LO[j])
    
    elig_below <- colSums(art_elig[j,,,drop=FALSE],,2) * pr_below
    if(j < hDS) 
      elig_below <- elig_below + colSums(art_elig[(j+1):hDS,,,drop=FALSE],,2)
    
    elig_above <- colSums(art_elig[j,,,drop=FALSE],,2) * (1.0-pr_below)
    if(j > 1) 
      elig_above <- elig_above + colSums(art_elig[1:(j-1),,,drop=FALSE],,2)

    initpr_below <- pmin(art_need * 0.5 / elig_below, 1.0, na.rm=TRUE)
    initpr_above <- pmin(art_need * 0.5 / elig_above, 1.0, na.rm=TRUE)
    initpr_medcat <- initpr_below * pr_below + initpr_above * (1-pr_below)

    out <- array(0, dim=c(hDS, hAG, NG))

    if(j < hDS)
      out[(j+1):hDS,,] <- sweep(art_elig[(j+1):hDS,,,drop=F], 3, initpr_below, "*")
    
    out[j,,] <- sweep(art_elig[j,,,drop=F], 3, initpr_medcat, "*")
    
    if(j > 0)
      out[1:(j-1),,] <- sweep(art_elig[1:(j-1),,,drop=F], 3, initpr_above, "*")
  }
  return(out)
}
# f_artDist(fp, art15plus.elig, art15plus.inits, i)

# Direct incidence input
f_infections_directincid <- function(pop, year, fp) {
  list2env(fp$ss, environment())
  if (fp$incidpopage == 0L) # incidence for 15 -49 population
    age_id <- p.age15to49.idx
  else if (fp$incidpopage == 1L) # incidence for 15+ population
    age_id <- p.age15plus.idx
  num <- c(1, fp$incrr_sex[year]) * sum(pop$data[age_id,, hivn.idx, year-1])
  den <- sum(pop$data[age_id, m.idx, hivn.idx, year-1]) + 
         sum(pop$data[age_id, f.idx, hivn.idx, year-1]) * fp$incrr_sex[year]
  sexinc <- fp$incidinput[year] * num /den
  ageinc <- colSums(pop$data[age_id,,hivn.idx,year-1] * fp$incrr_age[age_id,,year]) /
            colSums(pop$data[age_id,,hivn.idx,year-1])
  out    <- sweep(fp$incrr_age[,,year], 2, sexinc / ageinc, "*")
  return(out)
}

# Inside calc_infections_eppspectrum
# -----------------------------------------------------------------------------
f_art.ii <- function(HIV, ACTIVE, h, p, DT, i, ii, pop, hivpop, artpop, fp) {

  wDT <- (1-DT*(ii-1))
  art.ii <- sum(artpop[,,h,,i])

  if (sum(hivpop[,h[1],,i]) + sum(artpop[,,h[1],,i])  > 0) {
    art.ii <- art.ii - sum(pop[p[1], HIV & ACTIVE, i] * colSums(artpop[,,h[1],,i],,2) / (colSums(hivpop[,h[1],,i],,1) + colSums(artpop[,,h[1],,i],,2))) * wDT
  }

  if (sum(hivpop[,tail(h, 1)+1,,i]) + sum(artpop[,,tail(h, 1)+1,,i]) > 0) {
    art.ii <- art.ii + sum(pop[tail(p,1)+1, HIV & ACTIVE, i] * colSums(artpop[,,tail(h, 1)+1,,i],,2) / (colSums(hivpop[,tail(h, 1)+1,,i],,1) + colSums(artpop[,,tail(h, 1)+1,,i],,2))) * wDT
  }
  return(art.ii)
}

f_agesex.inc <- function(HIV, WOMAN, ACTIVE, inc.rate, p, i, pop, fp) {

  sexinc15to49.ts <- inc.rate * c(1, fp$incrr_sex[i]) * sum(pop[p, !HIV & ACTIVE, i]) / 
    (sum(pop[p, !WOMAN & !HIV & ACTIVE, i]) + fp$incrr_sex[i] * sum(pop[p, WOMAN & !HIV & ACTIVE, i]))

  sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts / (colSums(pop[p, !HIV & ACTIVE, i] * fp$incrr_age[p,,i]) / colSums(pop[p, !HIV & ACTIVE, i])), "*")
}

f_hiv.ii <- function(ages, HIV, ACTIVE, pop, fp, i, ii,  DT) {
  out <- sum(pop[ages, HIV & ACTIVE, i])
  out <- out - sum(pop[ages[1], HIV & ACTIVE, i]) * (1 - DT*(ii-1))
  out <- out + sum(pop[tail(ages,1)+1, HIV & ACTIVE, i]) * (1 - DT*(ii-1))
  return(out)
}

# End bits from eppasm.R
# Continue in progression.R
# -----------------------------------------------------------------------------