# Update HIV state space for debut model
update_hiv_ss <- function(max_debut_age, fp) {
  tt        <- fp$ss

  # recode the age indices for 
  AGE_END   <- tt$AGE_START + tt$pAG - 1
  tt$db_agr <- 1:tt$pDB
  agr_up    <- tt$ag.idx[-(tt$db_agr)]
  agr_up    <- tail(tt$db_agr, 1) + (agr_up - agr_up[1] + 1)
  agegr     <- c(tt$db_agr, agr_up)

  tt$ag.idx_bk       <- tt$ag.idx
  tt$ag.idx          <- as.integer(agegr)
  tt$h.ag.span       <- as.integer(rle(agegr)$lengths)
  tt$hAG             <- length(unique(agegr))
  tt$agfirst.idx     <- as.integer(which(!duplicated(tt$ag.idx)))
  tt$aglast.idx      <- as.integer(which(!duplicated(tt$ag.idx, fromLast=TRUE)))
  tt$h.fert.idx_bk   <- tt$h.fert.idx
  tt$h.fert.idx      <- as.integer(which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) %in% 15:49))
  tt$h.age15to49.idx <- as.integer(which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) %in% 15:49))
  tt$h.age15plus.idx <- as.integer(which((tt$AGE_START-1 + cumsum(tt$h.ag.span)) >= 15))

  # Expanding fp for the newly added age for sexual debut model
  fp$cd4_initdist <- fp$cd4_initdist[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$cd4_prog     <- fp$cd4_prog[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$cd4_mort     <- fp$cd4_mort[, tt$ag.idx_bk, ][, tt$agfirst.idx, ]
  fp$art_mort     <- fp$art_mort[,,tt$ag.idx_bk, ][,,tt$agfirst.idx, ]
  
  fertid_new   <- tt$ag.idx[tt$ag.idx %in% tt$h.fert.idx]
  agfirst_fert <- which(!duplicated(fertid_new))
  fertid       <- tt$ag.idx_bk[tt$ag.idx_bk %in% tt$h.fert.idx_bk]
  fp$frr_cd4   <- fp$frr_cd4[, fertid,][, agfirst_fert,]
  fp$frr_art   <- fp$frr_art[,,fertid,][,,agfirst_fert,]

  fp$ss <- tt
  return(fp)
}

#' Update fp for sexual debut transistion
#' 
#' @param fp fix parameters
#' @param max_debut_age all will become sexual active at this age
update_fp_debut <- function(fp, max_debut_age = 30) {
  if (!exists("db_rate", where=fp)) {
    cat('running with default sexual debut rate...\n')
    fp$db_rate     <- cbind(db_rate(fp$ss$AGE_START:max_debut_age),
                          db_rate(fp$ss$AGE_START:max_debut_age, FALSE))    
  }
  fp$db_rate[a2i(max_debut_age),,] <- 1
  fp$ss$db_aid <- a2i(fp$ss$AGE_START):a2i(max_debut_age)
  fp$ss$pDB    <- max(fp$ss$db_aid)
  fp           <- update_hiv_ss(max_debut_age, fp)
  return(fp)
}

#' Calculate debut rate for age and sex using log-logistic distribution
#' 
#' Note that the hazard is parameterize as in r-inla
#' 
#' @param age current age
#' @param male are we calculating for male or female?
#' @param shape the shape parameter for log-logistic distribution
#' @param scale the scale parameter for log-logistic distribution
#' @return Hazard rate
db_rate <- function(age, male=TRUE, shape, scale) {
  # average estimates from INLA for all countries
  if (missing(shape)) shape <- 8.849903
  if (missing(scale)) scale <- ifelse(male, 0.05338877, 0.05924315)
  hz_llogis(age, shape, scale)
}