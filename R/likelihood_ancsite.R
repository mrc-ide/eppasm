
#' Sample from posterior ditribution for ANC site level random effects
#'
sample_b_site <- function(mod, fp, dat, resid=TRUE){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)
  vinfl <- fp$v.infl

  df <- dat$df
  qM <- suppressWarnings(stats::qnorm(agepregprev(mod, fp, dat$datgrp$aidx, dat$datgrp$yidx, dat$datgrp$agspan)))

  mu <- qM[df$qMidx] + dat$Xancsite %*% coef
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
  qM <- suppressWarnings(stats::qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- stats::model.matrix(~type, df)

  mu <- qM[df$qMidx] + Xancsite %*% coef + b_site[match(df$site, names(b_site))]
  v <- 2 * pi * exp(mu^2) * stats::pnorm(mu) * (1 - stats::pnorm(mu))/df$n + fp$v.infl
  stats::rnorm(nrow(df), mu, sqrt(v))
}


#' Pointwise likelihood for site-level ANC observations given site-level effects
#'
ll_ancsite_conditional <- function(mod, fp, newdata, b_site){

  coef <- c(fp$ancbias, fp$ancrtsite.beta)

  df <- newdata$df
  qM <- suppressWarnings(stats::qnorm(agepregprev(mod, fp, newdata$datgrp$aidx, newdata$datgrp$yidx, newdata$datgrp$agspan)))

  ## Design matrix for fixed effects portion
  df$type <- factor(df$type, c("ancss", "ancrt"))
  Xancsite <- stats::model.matrix(~type, df)

  mu <- qM[df$qMidx] + Xancsite %*% coef + b_site[match(df$site, names(b_site))]
  v <- 2 * pi * exp(mu^2) * stats::pnorm(mu) * (1 - stats::pnorm(mu))/df$n + fp$v.infl
  stats::dnorm(df$W, mu, sqrt(v), log=TRUE)
}
