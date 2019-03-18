expand_age <- function(x) {
  apply(x[, 1:2], 2, rep, times = apply(x[, 3:4], 1, diff) + 1)
}

update_hiv_ss <- function(debut_pr, fp) {
  tt <- fp$ss
  AGE_END <- tt$AGE_START + tt$pAG - 1
  tt$db_agr <- unname(sapply(debut_pr[, 3], a2i))
  agr_up <- tt$ag.idx[-(tt$db_agr)]
  agr_up <- tail(tt$db_agr, 1) + (agr_up - agr_up[1] + 1)
  agegr  <- c(tt$db_agr, agr_up)

  tt$ag.idx_bk <- tt$ag.idx
  tt$ag.idx <- agegr
  tt$h.ag.span <- rle(agegr)$lengths
  tt$hAG <- length(unique(agegr))
  tt$agfirst.idx <- which(!duplicated(tt$ag.idx))
  tt$aglast.idx <- which(!duplicated(tt$ag.idx, fromLast=TRUE))
  tt$h.fert.idx_bk <- tt$h.fert.idx
  tt$h.fert.idx <- which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) %in% 15:49)
  tt$h.age15to49.idx <- which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) %in% 15:49)
  tt$h.age15plus.idx <- which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) >= 15)

  # Expanding fp
  fp$cd4_initdist <- fp$cd4_initdist[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$cd4_prog <- fp$cd4_prog[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$cd4_mort <- fp$cd4_mort[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$art_mort <- fp$art_mort[,,tt$ag.idx_bk, ][,,tt$agfirst.idx, ]
  
  fertid_new <- tt$ag.idx[tt$ag.idx %in% tt$h.fert.idx]
  agfirst_fert <- which(!duplicated(fertid_new))
  fertid <- tt$ag.idx_bk[tt$ag.idx_bk %in% tt$h.fert.idx_bk]
  fp$frr_cd4  <- fp$frr_cd4[, fertid,][, agfirst_fert,]
  fp$frr_art  <- fp$frr_art[,,fertid,][,,agfirst_fert,]

  fp$ss <- tt
  return(fp)
}

#' Update fp for sexual debut transistion
#' 
#' @param fp fix parameters
#' @param debut_pr data frame contains probability of sexual debut
#' @param rlq_pr data frame contains probability of sexual relinquish
#' @param single_year is the probability provided as single-year or age-group?
update_fp_debut <- function(fp, debut_pr, rlq_pr=NULL, single_year=FALSE) {
  list2env(fp$ss, environment())
  if (single_year) {
    # check for correct inputs
    if (!all.equal(debut_pr[,3], debut_pr[,4])) stop("Age male != age female")
    if (any(diff(debut_pr[,3]) > 1) | any(diff(debut_pr[,4]) > 1)) 
      stop("Single age is not continous")
    debut_pr <- debut_pr[debut_pr[,3] >= AGE_START,]
    ages <- debut_pr[, 3]
    fp$db_pr <- as.matrix(debut_pr[, 1:2])
    fp$ss$db_aid <- a2i(min(ages), AGE_START):a2i(max(ages), AGE_START)
    fp$ss$pDB <- max(fp$ss$db_aid)
    fp <- update_hiv_ss(debut_pr, fp)
    
    # Relinquish: this does not matter now as the prev is only to 15-49 but the
    # index need to be provided for when changing the prev
    if (is.null(rlq_pr)) {
      fp$rlq_pr   <- matrix(numeric(2), ncol=2)
      fp$ss$rlq_aid <- a2i(80, fp$ss$AGE_START)
    } else {
      fp$db_pr   <- rlq_pr[, 1:2]
      fp$ss$db_aid <- min(rlq_pr[, 3]):80 %>% sapply(a2i, fp$ss$AGE_START)
      fp$ss$pDB <- max(fp$ss$db_aid)
    }

  } else {
    debut_pr <- debut_pr[debut_pr[,3] >= AGE_START,]
    fp$db_pr <- expand_age(debut_pr)
    ages <- range(debut_pr[, 3:4])
    fp$ss$db_aid <- a2i(min(ages), AGE_START):a2i(max(ages), AGE_START)
    fp$ss$pDB <- max(fp$ss$db_aid)
    debut_pr <- cbind(fp$db_pr, ages[1]:ages[2], ages[1]:ages[2])
    fp <- update_hiv_ss(debut_pr, fp)

    # Relinquish: this does not matter now as the prev is only to 15-49
    if (is.null(rlq_pr)) {
      fp$rlq_pr   <- matrix(numeric(2), ncol=2)
      fp$ss$rlq_aid <- a2i(80, fp$ss$AGE_START)
    } else {
      fp$rlq_pr <- expand_age(rlq_pr)
      fp$ss$rlq_aid <- min(rlq_pr["age_lo"]):80 %>% sapply(a2i, fp$ss$AGE_START)
    }
  }
  return(fp)
}

# Put here for references :::
debut_pr <- data.frame(male = c(0.01, 0.27, 0.82, 1.00),
                     female = c(0.05, 0.52, 0.91, 1.00),
                     age_lo = c(10, 15, 20, 25),
                     age_up = c(14, 19, 24, 29))

# Relinquish: this does not matter now as the prev is only to 15-49
rlq_pr <- data.frame(male = c(.5, .8, 1),
                     female = c(.6, .9, 1),
                     age_lo = c(50, 55, 60),
                     age_up = c(54, 59, 80))

#' Distributing basepop
#' 
#' @param HIV logical vector generated from selectmodel()
#' @param ACTIVE logical vector generated from selectmodel()
#' @param SEX logical vector generated from selectmodel()
#' @param pop Core population
#' @param fp fix parameters
distribute_debut <- function(HIV, ACTIVE, SEX, pop, fp) {
  list2env(fp$ss, environment())  
  # Put to active 
  pop[, !HIV & ACTIVE, 1] <- fp$basepop
  # then move
  pop[db_aid, !HIV & !ACTIVE, 1] <- pop[db_aid, !HIV & ACTIVE, 1] * (1 - fp$db_pr)
  # the remain
  pop[db_aid, !HIV & ACTIVE, 1] <- pop[db_aid, !HIV & ACTIVE, 1] * fp$db_pr

  # move the olds
  pop[rlq_aid, !ACTIVE, 1] <- pop[rlq_aid, ACTIVE, 1] * fp$rlq_pr[, SEX[ACTIVE]+1 ]
  pop[rlq_aid,  ACTIVE, 1] <- pop[rlq_aid,  ACTIVE, 1] * (1 - fp$rlq_pr[, SEX[ACTIVE]+1]) # the remain
  return(pop)
}

calc_entrant <- function(MODEL, birth_lag, preg_prev_lag, fp, i) {
  if(exists("popadjust", where = fp) & fp$popadjust) {
    if (MODEL == 0) 
      return(fp$entrantpop[,i-1])
    if (MODEL != 0) {
      hivn <- fp$entrantpop[,i-1] * (1 - fp$entrantprev[,i])
      hivp <- fp$entrantpop[,i-1] * fp$entrantprev[,i]
      return(list(hivn, hivp))
    }
  } else {
    if (MODEL == 0) 
      return(HpIn0(fp, i, birth_lag))
    if (MODEL != 0) {
      hivn <- HnIn(fp, i, preg_prev_lag, birth_lag)
      hivp <- HpIn(fp, i, preg_prev_lag, birth_lag)
      return(list(hivn, hivp))
    }
  }
}

entrant_hiv <- function(HIV, ACTIVE, pop, hiv_pop, aging_prob, year, fix_par) {
  out <- age_hivpop(aging_prob[, HIV & ACTIVE], hiv_pop, year, fix_par)
  noART <- pop[1, HIV & ACTIVE, year] * (1-fix_par$entrantartcov[,year])
  out[,1,] <- out[,1,] + sweep(fix_par$paedsurv_cd4dist[,,year], 2, noART, "*")
  return(out)
}

entrant_art <- function(HIV, ACTIVE, pop, art_pop, aging_prob, year, fix_par) {
  out <- age_artpop(aging_prob[, HIV & ACTIVE], art_pop, year, fix_par)
  isART  <- pop[1, HIV & ACTIVE, year] * fix_par$entrantartcov[,year] # entrants treatment+
  out[,,1,] <- out[,,1,] + sweep(fix_par$paedsurv_artcd4dist[,,,year], 3, isART, "*")
  return(out)
}