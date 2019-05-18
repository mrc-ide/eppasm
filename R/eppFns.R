# -----------------------------------------------------------------------------
# Extracting bits from eppasm.R
# -----------------------------------------------------------------------------

# pop x by age groups
# When x is a vector
#' @importFrom fastmatch ctapply
sumByAG <- function(x, ag.idx, fertile=FALSE, p.fert.idx=NULL) {
  if (!fertile) fastmatch::ctapply(x, ag.idx, sum) 
    else fastmatch::ctapply(x, ag.idx[p.fert.idx], sum)
}

# pop x by age groups
# When x is a data.frame/matrix
sumByAGs <- function(k, ag.idx, fertile=FALSE, p.fert.idx=NULL) 
  apply(k, 2, sumByAG, ag.idx=ag.idx, fertile=fertile, p.fert.idx=p.fert.idx)

# Scale mortality cd4
scale_cd4_mort <- function(hivpop, artpop) {
  if (hivpop$p$scale_cd4_mort) {
    year <- hivpop$year
    num   <- hivpop$get(year) + hivpop$data_db[,,,year]
    den   <- colSums(artpop$get(year) + artpop$data_db[,,,,year])
    cd4mx <- num / (num + den)
    cd4mx[!is.finite(cd4mx)] <- 1.0
    cd4_mort_ts <- cd4mx * hivpop$p$cd4_mort
    return(cd4_mort_ts)
  } 
  else
    return(hivpop$p$cd4_mort)
}

# hiv deaths at ts
calc.agdist <- function(x, ag.idx, h.ag.span) {
  d <- x / rep(sumByAG(x, ag.idx), h.ag.span) # percent of each age/age group
  d[is.na(d)] <- 0 
  d
}

# End bits from eppasm.R
# Continue in progression.R
# -----------------------------------------------------------------------------