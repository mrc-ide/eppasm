#' Prepare ANC model object
#'
#' @details
#'
#' Defaults are to:
#' * Fit site-level ANC bias parameter if any ANC data exist
#' * Fit ANC-RT vs. ANC-SS site-level bias if ANC-RT site data exist
#' * Fit fertility rate ratio adjustment if ANC-RT census data exist
#' * Fit vinfl and ancrtcens_vinfl if ANC site or ANC data data exist, respectively
#' 
prepare_anc_model <- function(fp, eppd,
                              fit_vinfl = fp$fit_vinfl,
                              fit_ancrtcens_vinfl = fp$fit_ancrtcens_vinfl,
                              fit_ancrtsite_beta = fp$fit_ancrtsite_beta){

  ancmod <- list()

  ancmod$has_ancss <- !is.null(eppd$ancsitedat) && any(eppd$ancsitedat$type == "ancss")
  ancmod$has_ancrtsite <- !is.null(eppd$ancsitedat) && any(eppd$ancsitedat$type == "ancrt")
  ancmod$has_ancrtcens <- !is.null(eppd$ancrtcens) && nrow(eppd$ancrtcens)

  ancmod$fit_ancbias <- ancmod$has_ancss || ancmod$has_ancrt
  ancmod$fit_logfrr <- ancmod$has_ancrtcens
  ancmod$fit_vinfl <- if(is.null(fit_vinfl)) (ancmod$has_ancss || ancmod$has_ancrtsite) else fit_vinfl
  ancmod$fit_ancrtcens_vinfl <- if(is.null(fit_ancrtcens_vinfl)) ancmod$has_ancrtcens else fit_ancrtcens_vinfl
  ancmod$fit_ancrtsite_beta <- if(is.null(fit_ancrtsite_beta)) ancmod$has_ancss && ancmod$has_ancrtsite else fit_ancrtsite_beta
  
  if(!ancmod$has_ancrtsite)
    fp$ancrtsite.beta <- 0

  ancmod$nparam <- as.integer(ancmod$fit_ancbias + ancmod$fit_vinfl + ancmod$fit_logfrr + ancmod$fit_ancrtcens_vinfl + ancmod$fit_ancrtsite_beta)
  
  fp$ancmod <- ancmod
  fp
}

create_ancmod_param <- function(theta_ancmod, ancmod) {

  if(length(theta_ancmod) != ancmod$nparam)
    stop("Length of parameter vector not equal to expected number of parameters")

  param <- list()

  i <- 1
  if(ancmod$fit_ancbias){
    param$ancbias <- theta_ancmod[i]
    i <- i+1
   } else
     param$ancbias <- 0

  if(ancmod$fit_vinfl){
    param$v.infl <- exp(theta_ancmod[i])
    i <- i+1
   } else
     param$v.infl <- 0

  if(ancmod$fit_logfrr){
    param$log_frr_adjust <- theta_ancmod[i]
    i <- i+1
   } else
     param$log_frr_adjust <- 0

  if(ancmod$fit_ancrtcens_vinfl){
    param$ancrtcens.vinfl <- exp(theta_ancmod[i])
    i <- i+1
  } else
    param$ancrtcens.vinfl <- 0

  if(ancmod$fit_ancrtsite_beta){
    param$ancrtsite.beta <- theta_ancmod[i]
    i <- i+1
  } else
    param$ancrtsite.beta <- 0

  param
}

lprior_ancmod <- function(theta_ancmod, ancmod, prior_args) {

  for(i in seq_along(prior_args))
    assign(names(prior_args)[i], prior_args[[i]])

  lpr <- 0

  i <- 1
  if(ancmod$fit_ancbias){
    lpr <- lpr + dnorm(theta_ancmod[i], ancbias.pr.mean, ancbias.pr.sd, log=TRUE)
    i <- i+1
  }
  if(ancmod$fit_vinfl){
    lpr <- lpr + dexp(exp(theta_ancmod[i]), vinfl.prior.rate, TRUE) + theta_ancmod[i]
    i <- i+1
  }
  if(ancmod$fit_logfrr){
    lpr <- lpr + dnorm(theta_ancmod[i], log_frr_adjust.pr.mean, log_frr_adjust.pr.sd, log=TRUE)
    i <- i+1
  }
  if(ancmod$fit_ancrtcens_vinfl){
    lpr <- lpr + dexp(exp(theta_ancmod[i]), ancrtcens.vinfl.pr.rate, TRUE) + theta_ancmod[i]
    i <- i+1
  }
  if(ancmod$fit_ancrtsite_beta){
    lpr <- lpr + dnorm(theta_ancmod[i], ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd, log=TRUE)
    i <- i+1
  } 
  
  lpr
}

sample_prior_ancmod <- function(n, ancmod, prior_args) {

  for(i in seq_along(prior_args))
    assign(names(prior_args)[i], prior_args[[i]])
  
  mat <- matrix(NA, n, ancmod$nparam)

  i <- 1
  if(ancmod$fit_ancbias){
    mat[ , i] <- rnorm(n, ancbias.pr.mean, ancbias.pr.sd)
    i <- i+1
  }
  if(ancmod$fit_vinfl){
    mat[ , i] <- log(rexp(n, vinfl.prior.rate))
    i <- i+1
  }
  if(ancmod$fit_logfrr){
    mat[ , i] <- rnorm(n, log_frr_adjust.pr.mean, log_frr_adjust.pr.sd)
    i <- i+1
  }
  if(ancmod$fit_ancrtcens_vinfl){
    mat[ , i] <- log(rexp(n, ancrtcens.vinfl.pr.rate))
    i <- i+1
  }
  if(ancmod$fit_ancrtsite_beta){
    mat[ , i] <- rnorm(n, ancrtsite.beta.pr.mean, ancrtsite.beta.pr.sd)
    i <- i+1
  } 

  mat
}


###########################################
####                                   ####
####  Site level ANC data (SS and RT)  ####
####                                   ####
###########################################

ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0
vinfl.prior.rate <- 1/0.015

ancrtsite.beta.pr.mean <- 0
ancrtsite.beta.pr.sd <- 1.0


#' Prepare design matrix indices for ANC prevalence predictions
#'
#' @param ancsite_df data.frame of site-level ANC design for predictions
#' @param fp fixed parameter input list
#'
#' @examples
#' pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
#' bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
#'
#' 
#' bw_u_ancsite <- attr(bw$Urban, "eppd")$ancsitedat
#' fp <- attr(bw$Urban, "specfp")
#'
#' ancsite_pred_df(bw_u_ancsite, fp)
#' 
ancsite_pred_df <- function(ancsite_df, fp) {

  df <- ancsite_df
  anchor.year <- fp$ss$proj_start
  
  df$aidx <- df$age - fp$ss$AGE_START + 1L
  df$yidx <- df$year - anchor.year + 1
  
  ## List of all unique agegroup / year combinations for which prevalence is needed
  datgrp <- unique(df[c("aidx", "yidx", "agspan")])
  datgrp$qMidx <- seq_len(nrow(datgrp))

  ## Indices for accessing prevalence offset from datgrp
  df <- merge(df, datgrp)
  
  list(df = df, datgrp = datgrp)
}


#' Prepare site-level ANC prevalence data for EPP random-effects likelihood
#'
#' @param ancsitedat data.frame of site-level ANC data
#' @param fp fixed parameter input list, including state space
#' @param offset column name for probit ANC prevalence offset 

prepare_ancsite_likdat <- function(ancsitedat, fp, offset = "offset"){

  d <- ancsite_pred_df(ancsitedat, fp)

  df <- d$df
  
  ## Calculate probit transformed prevalence and variance approximation
  df$pstar <- (df$prev * df$n + 0.5) / (df$n + 1)
  df$W <- qnorm(df$pstar)
  df$v <- 2 * pi * exp(df$W^2) * df$pstar * (1 - df$pstar) / df$n

  ## offset
  if(length(offset) > 1 ||
     !is.character(offset) ||
     offset != "offset" && is.null(df[[offset]]))
    stop("ANC offset column not found")
  else if(is.null(df[[offset]]))
    df$offset <- 0
  else 
    df$offset <- df[[offset]]
  
  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- model.matrix(~type, df)

  ## Indices for observation
  df_idx.lst <- split(seq_len(nrow(df)), factor(df$site))

  df <- df[c("site", "year", "used", "type", "age", "agspan",
             "n", "prev", "aidx", "yidx", "qMidx",
             "pstar", "W", "v", "offset", "type")]

  list(df = df,
       datgrp = d$datgrp,
       Xancsite = Xancsite,
       df_idx.lst = df_idx.lst)
}

ll_ancsite <- function(mod, fp, coef=c(0, 0), vinfl=0, dat){

  df <- dat$df

  if(!nrow(df))
    return(0)

  if(exists("pregprev", fp) && !fp$pregprev)
    qM <- suppressWarnings(qnorm(ageprev(mod, dat$datgrp$aidx, rep(0L, nrow(dat$datgrp)), dat$datgrp$yidx, dat$datgrp$agspan)))
  else
    qM <- suppressWarnings(qnorm(agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx, dat$datgrp$agspan)))

  if(any(is.na(qM)) || any(qM == -Inf) || any(qM > 2))  ## prev < 0.977
    return(-Inf)
  
  mu <- qM[df$qMidx] + dat$Xancsite %*% coef + df$offset
  d <- df$W - mu
  v <- df$v + vinfl
  
  d.lst <- lapply(dat$df_idx.lst, function(idx) d[idx])
  v.lst <- lapply(dat$df_idx.lst, function(idx) v[idx])
  
  log(anclik::anc_resid_lik(d.lst, v.lst))
}


#' Sample from posterior ditribution for ANC site level random effects
#' 
sample_b_site <- function(mod, fp, dat, resid=TRUE){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)
  vinfl <- fp$v.infl

  df <- dat$df

  if(exists("pregprev", fp) && !fp$pregprev)
    qM <- suppressWarnings(qnorm(ageprev(mod, newdata$aidx, rep(0L, nrow(newdata)), newdata$yidx, newdata$agspan)))
  else
    qM <- suppressWarnings(qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  mu <- qM[df$qMidx] + dat$Xancsite %*% coef + df$offset
  d <- df$W - mu
  v <- df$v + vinfl
  
  d.lst <- lapply(dat$df_idx.lst, function(idx) d[idx])
  v.lst <- lapply(dat$df_idx.lst, function(idx) v[idx])

  b_site <- mapply(anclik::sample.b.one, d.lst, v.lst)

  if(resid){
    mu_site <- mu + b_site[match(df$site, names(b_site))]
    resid <- c(df$W - mu_site)
    return(list(b_site = b_site, resid = resid))
  }

  b_site
}


#' Sample posterior predictions for site-level ANC observations
#' 
sample_ancsite_pred <- function(mod, fp, newdata, b_site){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)
  
  df <- newdata$df

  if(exists("pregprev", fp) && !fp$pregprev)
    qM <- suppressWarnings(qnorm(ageprev(mod, newdata$aidx, rep(0L, nrow(newdata)), newdata$yidx, newdata$agspan)))
  else
    qM <- suppressWarnings(qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- model.matrix(~type, df)

  mu <- qM[df$qMidx] + Xancsite %*% coef + df$offset + b_site[match(df$site, names(b_site))]
  v <- 2 * pi * exp(mu^2) * pnorm(mu) * (1 - pnorm(mu))/df$n + fp$v.infl
  rnorm(nrow(df), mu, sqrt(v))
}  


#' Pointwise likelihood for site-level ANC observations given site-level effects
#' 
ll_ancsite_conditional <- function(mod, fp, newdata, b_site){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)
  
  df <- newdata$df

  if(exists("pregprev", fp) && !fp$pregprev)
    qM <- suppressWarnings(qnorm(ageprev(mod, newdata$aidx, rep(0L, nrow(newdata)), newdata$yidx, newdata$agspan)))
  else
    qM <- suppressWarnings(qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- model.matrix(~type, df)

  mu <- qM[df$qMidx] + Xancsite %*% coef + df$offset + b_site[match(df$site, names(b_site))]
  v <- 2 * pi * exp(mu^2) * pnorm(mu) * (1 - pnorm(mu))/df$n + fp$v.infl
  dnorm(df$W, mu, sqrt(v), log=TRUE)
}  


#############################################
####                                     ####
####  ANCRT census likelihood functions  ####
####                                     ####
#############################################

## prior parameters for ANCRT census
log_frr_adjust.pr.mean <- 0
log_frr_adjust.pr.sd <- 1.0
ancrtcens.vinfl.pr.rate <- 1/0.015

prepare_ancrtcens_likdat <- function(dat, fp){

  anchor.year <- fp$ss$proj_start
  
  x.ancrt <- (dat$prev*dat$n+0.5)/(dat$n+1)
  dat$W.ancrt <- qnorm(x.ancrt)
  dat$v.ancrt <- 2*pi*exp(dat$W.ancrt^2)*x.ancrt*(1-x.ancrt)/dat$n

  if(!exists("age", dat))
    dat$age <- rep(15, nrow(dat))

  if(!exists("agspan", dat))
    dat$agspan <- rep(35, nrow(dat))

  dat$aidx <- dat$age - fp$ss$AGE_START + 1
  dat$yidx <- dat$year - anchor.year + 1
  
  return(dat)
}

ll_ancrtcens <- function(mod, dat, fp, pointwise = FALSE){
  if(!nrow(dat))
    return(0)

  if(exists("pregprev", fp) && !fp$pregprev)
    qM <- suppressWarnings(qnorm(ageprev(mod, dat$aidx, rep(0L, nrow(dat)), dat$yidx, dat$agspan)))
  else
    qM <- suppressWarnings(qnorm(agepregprev(mod, fp, dat$aidx, dat$yidx, dat$agspan)))

  if(any(is.na(qM)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- dnorm(dat$W.ancrt, qM, sqrt(dat$v.ancrt + fp$ancrtcens.vinfl), log=TRUE)

  if(pointwise)
    return(val)

  sum(val)
}
