outpred_prev <- function(fit, out, subsample=NULL){

  hhs <- out$hhs

  if(!is.null(subsample))
    fit$resample <- fit$resample[sample.int(nrow(fit$resample), subsample), ]

  dat <- prepare_hhsageprev_likdat(hhs, fit$fp)
  fit <- simfit(fit, ageprevdat=dat, pregprev=FALSE, ageincid=FALSE, ageinfections=FALSE,
                relincid=FALSE, entrantprev=FALSE)

  ## Impute missing standard errors based on model prevalence
  na_i <- which(is.na(dat$sd.W.hhs))
  na_p <- rowMeans(fit$ageprevdat)[na_i]
  dat$sd.W.hhs[na_i] <- sqrt(na_p * (1 - na_p) / dat$n_eff[na_i]) / dnorm(qnorm(na_p))


  M <- fit$ageprevdat
  qM <- qnorm(M)
  qpred <- array(rnorm(length(qM), qM, dat$sd.W.hhs), dim(qM))

  elpd <- log(rowMeans(exp(ldbinom(dat$x_eff, dat$n_eff, M))))
  resid <- rowMeans(M) - dat$prev
  rse <- apply(M, 1, sd) / rowMeans(M)
  mae <- rowMeans(abs(M - dat$prev))
  rmse <- sqrt(rowMeans((M - dat$prev)^2))

  elpd_q <- log(rowMeans(dnorm(dat$W.hhs, qM, dat$sd.W.hhs)))
  rse_q <- apply(qM, 1, sd) / rowMeans(qM)
  resid_q <- rowMeans(qM) - dat$W.hhs
  mae_q <- rowMeans(abs(qM - dat$W.hhs))
  rmse_q<- sqrt(rowMeans((qM - dat$W.hhs)^2))
  qq <- mapply(function(f, x) f(x), apply(qpred, 1, ecdf), dat$W.hhs)  

  vars <- intersect(c("country", "eppregion", "survyear", "year", "sex", "agegr", "prev", "se"), names(dat))
  data.frame(dat[vars],
             elpd, resid, rse, mae, rmse, elpd_q, resid_q, rse_q, mae_q, rmse_q, qq)
}
