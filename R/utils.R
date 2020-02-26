sweepx <- function(...) sweep(..., FUN='*') # missing FUN too many times!

# Natural age to index
a2i <- function(x, min=15, max=80) which(min:max %in% x)

# Non number to value
na2num <- function(x, y) {x[is.na(x)] <- y; return(x)}

# hazard function for log-logistic distribution (parameterize as in INLA)
hz_llogis <- function (x, alpha = 8, lambda = 1/18) {
    num <- alpha * lambda
    den <- (lambda * x)^(1-alpha) + lambda * x
    num/den
}
# Cumulative function for log-logistic distribution (parameterize as in INLA)
cdf_llogis <- function (x, alpha = 8, lambda = 1/18) {
    1 / ( 1 + (lambda * x)^(-alpha) )
}

logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))
ldinvlogit <- function(x){v <- invlogit(x); log(v) + log(1-v)}

#' Replace R naming of fp for C++ 
#' 
#' dots are replaced with underscore, char is replaced with integer, 
#' set default value for mixing model and dummy variables
#' 
#' @param fp Fix parameters
prepare_fp_for_Cpp <- function(fp, MODEL=1L, MIX=FALSE) {
    # This is a mess
    # rhybrid = 0 # rtrend = 1 # directincid =2
    fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid"), nomatch=0)
    fp$ss$MODEL <- as.integer(MODEL)
    fp$ss$MIX   <- MIX
    names(fp$ss) <- gsub('\\.', '_', names(fp$ss))
    names(fp) <- gsub('\\.', '_', names(fp))
    if (exists("rt", where=fp))
        names(fp$rt) <- gsub('\\.', '_', names(fp$rt))
    if (!exists("DT", where=fp))
        fp$ss$DT <- 1 / fp$ss$hiv_steps_per_year
    if ( is.integer(fp$ss$time_epi_start) )
        fp$ss$time_epi_start <- as.numeric(fp$ss$time_epi_start)
    fp$incidmodInt <- match(fp$incidmod, c("eppspectrum"))-1L
    fp$ancrtInt <- match(fp$ancrt, c("both"), nomatch=0) # just a placeholder
    fp$art15plus_isperc <- apply(fp$art15plus_isperc, 2, as.numeric)
    fp$ss$h_ag_span <- as.numeric(fp$ss$h_ag_span)
    if (!exists("rw_start", where=fp)) 
        fp$rw_start <- fp$rt$rw_start;
    if (!exists("rvec", where=fp)) 
        fp$rvec <- 1
    if (!exists("pDB", where=fp$ss))  { # mixing model parameters
        fp$ss$pDB <- 1L
        fp$db_rate <- array(1, c(fp$ss$pDB, 2, 52))
    }
    if (!exists("mixmat", where=fp))  { # mixing model parameters
        fp$mixmat <- array(0, c(fp$ss$pAG, fp$ss$pAG, fp$ss$NG))
    }
    if (!fp$popadjust) {
        fp$targetpop <- array(0, c(1,1,1))
        fp$entrantpop <- array(0, c(1,1))
    }
    if (is.null(dim(fp$artmx_timerr))) { # repicate for 3 treatment durations
        fp$artmx_timerr <- matrix(rep(fp$artmx_timerr, 3), nrow=3, byrow=TRUE)
    }
    if (fp$ss$MIX) { # scaled for sexual mixing model C++
        max_by_year <- apply(fp$incrr_age, 3, function(x) apply(x, 2, max))
        fp$incrr_age <- sweep(fp$incrr_age, 2:3, max_by_year, "/")
    }
    fp
}
# Converting prior assumption to parameter boundary for DE
prior_to_DE_bounds <- function(fp) {

  if (exists("prior_args", where = fp)){
	for(i in seq_along(fp$prior_args))
	  assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  up <- c()
  lo <- c()

  if (fp$eppmod == "rhybrid") {
	up <- rlog_pr_mean + 2.58*rlog_pr_sd
	lo <- rlog_pr_mean - 2.58*rlog_pr_sd
	up <- c(up, rep(+2.58*rw_prior_sd, fp$rt$n_rw))
	lo <- c(lo, rep(-2.58*rw_prior_sd, fp$rt$n_rw))
	if (exists("logitiota", fp) && fp$logitiota) {
		up <- c(up, logit(.9999))
		lo <- c(lo, logit(.0001))
	} else {
		up <- c(up, logiota.unif.prior[2])
		lo <- c(lo, logiota.unif.prior[1])
	}
  }
  ## sample ANC model parameters
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
	  if(fp$ancmod$fit_ancbias) {
		up <- c(up, ancbias.pr.mean + 2.58 * ancbias.pr.sd)
		lo <- c(lo, ancbias.pr.mean - 2.58 * ancbias.pr.sd)
	  }
	  if(fp$ancmod$fit_vinfl){
		bs <- range(log(rexp(1000, ancrtcens.vinfl.pr.rate)))
		lo <- c(lo, bs[1])
		up <- c(up, bs[2])
	  }
	  if(fp$ancmod$fit_logfrr){
		lo <- c(lo, log_frr_adjust.pr.mean - 2.58 * log_frr_adjust.pr.sd)
		up <- c(up, log_frr_adjust.pr.mean + 2.58 * log_frr_adjust.pr.sd)
	  }
	  if(fp$ancmod$fit_ancrtcens_vinfl){
		bs <- range(log(rexp(1000, ancrtcens.vinfl.pr.rate)))
		lo <- c(lo, bs[1])
		up <- c(up, bs[2])
	  }
	  if(fp$ancmod$fit_ancrtsite_beta){
	  	lo <- c(lo, ancrtsite.beta.pr.mean - 2.58 * ancrtsite.beta.pr.sd)
		up <- c(up, ancrtsite.beta.pr.mean + 2.58 * ancrtsite.beta.pr.sd)
	  } 
  }

  return(cbind(lo, up))
}
