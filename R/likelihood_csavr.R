

#' Basic logistic function for incidence rate
ilogistic <- function(t, p, t0){
  ## p[1] = alpha (growth rate)
  ## p[2] = c (max value)

  e <- exp(p[1] * (t - t0))
  p[2] * e / (1+e)
}

#' Double logistic function for incidence rate
idbllogistic <- function(t, p){
  ## p[1] = alpha
  ## p[2] = beta
  ## p[3] = t0
  ## p[4] = a
  ## p[5] = b
  e1 <- exp(p[1] * (t - p[3]))
  e2 <- exp(-p[2] * (t - p[3]))

  e1/(1+e1) * (2 * p[4] * (e2/(1+e2)) + p[5])
}


#' Calculate predicted number of new diagnoses
calc_diagnoses <- function(mod, fp){

  colSums(attr(mod, "diagnoses"),,3)
}

cumgamma_diagn_rate <- function(gamma_max, delta_rate, fp){
  
  delta_t <- rep(0, fp$ss$PROJ_YEARS)
  ii <- fp$t_diagn_start:fp$ss$PROJ_YEARS

  delta_t[ii] <- gamma_max * pgamma(ii - (ii[1] - 1), shape=1, rate = delta_rate)

  ## Diagnosis rate assumed proportional to expected mortality
  fp$diagn_rate <- array(fp$cd4_mort,
                         c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, fp$ss$PROJ_YEARS))

  ii <- seq_len(fp$ss$PROJ_YEARS)
  delta_t <- fp$gamma_max * pgamma(ii, shape=1, rate = fp$delta_rate)

  fp$diagn_rate <- sweep(fp$diagn_rate, 4, delta_t, "*")

  fp
}


#' Calculate parameter inputs for CSAVR fit
create_param_csavr <- function(theta, fp){
  
  if(fp$eppmod == "directincid" && fp$incid_func == "ilogistic"){
    nparam_incid <- 2
    fp$incidinput <- ilogistic(seq_len(fp$ss$PROJ_YEARS), exp(theta[1:2]), 1)
  }

  if(fp$eppmod == "directincid" && fp$incid_func == "idbllogistic"){
    nparam_incid <- 5
    tt <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
    p <- theta[1:5]
    p[c(1:2, 4:5)] <- exp(theta[c(1:2, 4:5)])
    fp$incidinput <- idbllogistic(tt, p)
  }

  if(fp$eppmod == "rlogistic"){
    nparam_incid <- 5
    par <- theta[1:4]
    par[3] <- exp(theta[3])
    fp$rvec <- exp(rlogistic(fp$proj.steps, par))
    fp$iota <- exp(theta[5])
  }

  if(fp$eppmod == "logrw"){
    nparam_incid <- fp$numKnots + 1L
    beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = exp(as.vector(fp$rvec.spldes %*% beta)),
                  iota = transf_iota(theta[fp$numKnots+1], fp))
    fp[names(param)] <- param
  }

  fp$gamma_max <- exp(theta[nparam_incid+1])
  fp$delta_rate <- exp(theta[nparam_incid+2])

  fp <- cumgamma_diagn_rate(fp$gamma_max, fp$delta_rate, fp)

  fp
}

#' Log-liklihood for new diagnoses and AIDS deaths
ll_csavr <- function(theta, fp, likdat){

  fp <- create_param_csavr(theta, fp)

  mod <- simmod(fp)

  mod_aidsdeaths <- colSums(attr(mod, "hivdeaths"),,2)
  ll_aidsdeaths <- with(likdat$aidsdeaths, sum(dpois(aidsdeaths, mod_aidsdeaths[idx] * (1 - prop_undercount), log=TRUE)))

  mod_diagnoses <- calc_diagnoses(mod, fp)
  ll_diagnoses <- with(likdat$diagnoses, sum(dpois(diagnoses, mod_diagnoses[idx] * (1 - prop_undercount), log=TRUE)))

  ll_aidsdeaths + ll_diagnoses
}



ilogistic_theta_mean <- c(-1, -10)
ilogistic_theta_sd <- c(5, 5)

idbllogistic_theta_mean <- c(-1, -1, 1995, -10, -10)
idbllogistic_theta_sd <- c(5, 5, 10, 5, 5)

diagn_theta_mean <- c(3, -3)
diagn_theta_sd <- c(5, 5)

logiota_pr_mean <- -13
logiota_pr_sd <- 5

sample_prior_csavr <- function(n, fp){

  mat_eppmod <- sample_prior_eppmod(n, fp)
  mat_diagn <- sample_prior_diagn(n, fp)

  cbind(mat_eppmod, mat_diagn)
}

sample_prior_eppmod <- function(n, fp){

  if(fp$eppmod == "logrw"){
    nparam <- fp$numKnots + 1L

    mat <- matrix(NA, n, nparam)
    mat[,1] <- rnorm(n, 0.2, 1)  # u[1]
    mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    mat[,fp$numKnots+1] <- sample_iota(n, fp)

  } else {

    theta_mean <- numeric()
    theta_sd <- numeric()

    if(fp$eppmod == "directincid" && fp$incid_func == "ilogistic"){
    theta_mean <- c(theta_mean, ilogistic_theta_mean)
    theta_sd <- c(theta_sd, ilogistic_theta_sd)
    }
    else if(fp$eppmod == "directincid" && fp$incid_func == "idbllogistic"){
      theta_mean <- c(theta_mean, idbllogistic_theta_mean)
      theta_sd <- c(theta_sd, idbllogistic_theta_sd)
    }
  else if(fp$eppmod == "rlogistic"){
      theta_mean <- c(theta_mean, rlog_pr_mean, logiota_pr_mean)
      theta_sd <- c(theta_sd, rlog_pr_sd, logiota_pr_sd)
    }
  else if(fp$eppmod == "rlogrw"){
      epp_nparam <- fp$numKnots+1L
    }
  else
    stop("incidence model not recognized")

  nparam <- length(theta_mean)

    ## Create matrix of samples
    v <- rnorm(nparam * n, theta_mean, theta_sd)
    mat <- t(matrix(v, nparam, n))
  }

  mat
}

sample_prior_diagn <- function(n, fp){

  nparam <- length(diagn_theta_mean)
  val <- rnorm(n * nparam, diagn_theta_mean, diagn_theta_sd)
  mat <- t(matrix(val, nparam, n))

  mat
}

lprior_csavr <- function(theta, fp){

  nparam_eppmod <- get_nparam_eppmod(fp)
  nparam_diagn <- 2L

  lprior_eppmod(theta[1:nparam_eppmod], fp) +
    lprior_diagn(theta[nparam_eppmod + 1:nparam_diagn])
}

lprior_eppmod <- function(theta_eppmod, fp){

  if(fp$eppmod == "directincid" && fp$incid_func == "ilogistic")
    return(sum(dnorm(theta_eppmod, ilogistic_theta_mean, ilogistic_theta_sd, log=TRUE)))
  else if(fp$eppmod == "directincid" && fp$incid_func == "idbllogistic")
    return(sum(dnorm(theta_eppmod, idbllogistic_theta_mean, idbllogistic_theta_sd, log=TRUE)))
  else if(fp$eppmod == "rlogistic")
    return(sum(dnorm(theta_eppmod, rlog_pr_mean, rlog_pr_sd, log=TRUE)))
  else if(fp$eppmod == "logrw"){
    lpr <- bayes_lmvt(theta_eppmod[2:fp$numKnots], rw_prior_shape, rw_prior_rate)
    lpr <- lpr + lprior_iota(theta_eppmod[fp$numKnots+1], fp)
    return(lpr)
  }
  else
    stop("incidence model not recognized")

}

lprior_diagn <- function(theta_diagn, fp){
  sum(dnorm(theta_diagn, diagn_theta_mean, diagn_theta_sd, log=TRUE))
}

get_nparam_eppmod <- function(fp){
  if(fp$eppmod == "directincid" && fp$incid_func == "ilogistic")
    return(length(ilogistic_theta_mean))
  else if(fp$eppmod == "directincid" && fp$incid_func == "idbllogistic")
    return(length(idbllogistic_theta_mean))
  else if(fp$eppmod == "rlogistic")
    return(length(rlog_pr_mean))
  else if(fp$eppmod == "logrw")
    return(fp$numKnots + 1L)
  else
    stop("incidence model not recognized")
}

prior_csavr <- function(theta, fp, log=FALSE){
  if(is.vector(theta))
    lval <- lprior_csavr(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (lprior_csavr(theta[i,], fp))))
  if(log)
    return(lval)
  else
    return(exp(lval))
}

likelihood_csavr <- function(theta, fp, likdat, log=FALSE){
  if(is.vector(theta))
    lval <- ll_csavr(theta, fp, likdat)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) ll_csavr(theta[i,], fp, likdat)))
  if(log)
    return(lval)
  else
    return(exp(lval))
}
