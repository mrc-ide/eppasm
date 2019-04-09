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

  ancmod$has_ancss <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancss")
  ancmod$has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
  ancmod$has_ancrtcens <- exists("ancrtcens", eppd) && nrow(eppd$ancrtcens)

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

#' Sample from posterior ditribution for ANC site level random effects
#' 
sample_b_site <- function(mod, fp, dat, resid=TRUE){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)
  vinfl <- fp$v.infl

  df <- dat$df
  qM <- suppressWarnings(qnorm(agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx, dat$datgrp$agspan)))

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
  qM <- suppressWarnings(qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- model.matrix(~type, df)

  mu <- qM[df$qMidx] + Xancsite %*% coef + df$offset + b_site[match(df$site, names(b_site))]
  v <- 2 * pi * exp(mu^2) * pnorm(mu) * (1 - pnorm(mu))/df$n + fp$v.infl
  dnorm(df$W, mu, sqrt(v), log=TRUE)
}  
