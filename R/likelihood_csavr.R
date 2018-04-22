

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

  ii <- seq_len(fp$ss$PROJ_YEARS)
  delta_t <- fp$gamma_max * pgamma(ii, shape=1, rate = fp$delta_rate)

  val <- attr(mod, "hivpop")
  val <- sweep(sweep(val, 1:3, fp$cd4_mort, "*"), 4, delta_t, "*")
  val <- colSums(val,, 3)

  return(val)
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
    
  fp$gamma_max <- exp(theta[nparam_incid+1])
  fp$delta_rate <- exp(theta[nparam_incid+2])

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
  else
    stop("incidnce model not recognized")

  theta_mean <- c(theta_mean, diagn_theta_mean)
  theta_sd <- c(theta_sd, diagn_theta_sd)

  nparam <- length(theta_mean)
  
  ## Create matrix of samples
  v <- rnorm(nparam * n, theta_mean, theta_sd)
  mat <- t(matrix(v, nparam, n))

  mat
}

lprior_csavr <- function(theta, fp){

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

  else
    stop("incidnce model not recognized")

  theta_mean <- c(theta_mean, diagn_theta_mean)
  theta_sd <- c(theta_sd, diagn_theta_sd)

  sum(dnorm(theta, theta_mean, theta_sd, log=TRUE))
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
