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
        fp$db_pr <- matrix(1, fp$ss$pDB, 2)
    }
    if (!exists("mat_f", where=fp))  { # mixing model parameters
        fp$mat_f <- fp$mat_m <- matrix(1, fp$ss$pAG, fp$ss$pAG)
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