

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
  
  if(fp$likelihood_cd4 == F){
  
    out <- colSums(attr(mod, "diagnoses"),,3)
    
  }
  
  if(fp$likelihood_cd4 == T){
  diagnoses_full <- attributes(mod)

  diagnoses <- diagnoses_full$diagnoses

  disease_stage <- as.character(1:nrow(diagnoses))

  diagnoses_df <- NULL
  for(i in 1:nrow(diagnoses)){
    sum_per_stage <- NULL
    for(ii in 1:52){
      sum_per_stage[ii] <- sum(diagnoses[i,,,ii])
    }
    sum_per_stage <- cbind.data.frame(sum_per_stage,c(1970:2021),rep(disease_stage[i],length(sum_per_stage)))
    diagnoses_df <- rbind.data.frame(diagnoses_df,sum_per_stage)
  }
  names(diagnoses_df) <- c("diagnoses","year","stage")
  
  stage_1_data <- diagnoses_df[diagnoses_df$stage == "1", ]
  names(stage_1_data) <- c("diagnoses","year","stage")
  
  stage_2_data <- diagnoses_df[diagnoses_df$stage == "2", ]
  names(stage_2_data) <- c("diagnoses","year","stage")
  
  stage_3_data <- cbind.data.frame((diagnoses_df[diagnoses_df$stage == "3", 1] + diagnoses_df[diagnoses_df$stage == "4", 1]),
                                   c(1970:2021), rep("3", nrow(stage_1_data)))
  names(stage_3_data) <- c("diagnoses","year","stage")
  
  stage_4_data <- cbind.data.frame((diagnoses_df[diagnoses_df$stage == "5", 1] + diagnoses_df[diagnoses_df$stage == "6", 1] + 
                                      diagnoses_df[diagnoses_df$stage == "7", 1]),
                                   c(1970:2021), rep("4", nrow(stage_1_data)))
  names(stage_4_data) <- c("diagnoses","year","stage")
  diagnoses_df_condensed <- rbind.data.frame(stage_1_data, stage_2_data, stage_3_data, stage_4_data)
  
  out <- diagnoses_df_condensed
  }

  return(out)
}

calc_artinits <- function(mod, fp){
  
  if(fp$likelihood_cd4 == F){
    
    art_df <- colSums(attr(mod,"artinits"),,3)
  }else{
  
  art_full <- attributes(mod)
  
  art <- art_full$artinits
  
  disease_stage_art <- as.character(1:nrow(art))
  
  art_df <- NULL
  for(i in 1:nrow(art)){
    tot_per_stage <- NULL
    for(ii in 1:52){
      tot_per_stage[ii] <- sum(art[i,,,ii])
    }
    
    tot_per_stage <- cbind.data.frame(tot_per_stage,c(1970:2021),rep(disease_stage_art[i],length(tot_per_stage)))
    art_df <- rbind.data.frame(art_df, tot_per_stage)
  }
  names(art_df) <- c("artinit","year","stage")
  
  stage_1_data <- art_df[art_df$stage == "1", ]
  names(stage_1_data) <- c("artinit","year","stage")
  
  stage_2_data <- art_df[art_df$stage == "2", ]
  names(stage_2_data) <- c("artinit","year","stage")
  
  stage_3_data <- cbind.data.frame((art_df[art_df$stage == "3", 1] + art_df[art_df$stage == "4", 1]),
                                c(1970:2021), rep("3", nrow(stage_1_data)))
  names(stage_3_data) <- c("artinit","year","stage")
  
  stage_4_data <- cbind.data.frame((art_df[art_df$stage == "5", 1] + art_df[art_df$stage == "6", 1] + 
                                      art_df[art_df$stage == "7", 1]),
                                   c(1970:2021), rep("4", nrow(stage_1_data)))
  names(stage_4_data) <- c("artinit","year","stage")
  
  art_df <- rbind.data.frame(stage_1_data, stage_2_data, stage_3_data, stage_4_data)
  }
  
  return(art_df)
  
}

knot_linear_diagn_rate <- function(knot_params, fp){
  
  # knots <- c(1970, 1980, 1986, 1996, 2001, 2009, 2015)
  # theta <- c(0, knot_params)
  knots <- c(1986, 1996, 2000, 2009, 2015)
  theta <- knot_params
  
  
  diagn_trend <- approx(knots, theta, 1970:2021, rule = 2)$y
  
  ## Mortality rate relative to Men, 25-34, CD4 100-200
  cd4_rel_diagn <- fp$cd4_mort / fp$cd4_mort[5, 4, 1]
  
  fp$diagn_rate <- array(cd4_rel_diagn,
                              c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, fp$ss$PROJ_YEARS))
  fp$diagn_rate <- sweep(fp$diagn_rate, 4, diagn_trend, "*")
  
  return(fp)
  
}

cumulative_linear_diagn_rate <- function(linear_params, fp){
  
  ## Linear params should be 12 values:
  ## First two represent increase rate for cd4 <350 from 1986 to 2000
  ## Second two represent baseline of cd4 > 350 from 1996 to 2000
  ## Params 5-8 represent linear increase for the detection rates of each cd4 class from 2001 to 2010 with SISCEL
  ## Params 9-12 represent linear increase for the detection rates of each cd4 class from 2011 to 2015 
  
  detec_array <- array(0, c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, fp$ss$PROJ_YEARS))
  
  ## base detection rate for the last two cd4 classes from 1981 to 1985
  
  detec_array[c(3:4),,,11:16] <- 0.01
  detec_array[5:7,,,11:16] <- 0.05
  
  ### linear increase 1986 - 2000
  
  detec_array[3:4,,,17:31] <- sweep(detec_array[3:4,,,17:31],4,(0.01 + seq(1,15,1) * linear_params[1]),"+")
  detec_array[5:7,,,17:31] <- sweep(detec_array[5:7,,,17:31],4,(0.05 + seq(1,15,1) * linear_params[2]),"+")
  
  ## For highest cd4 classes, detection in this instance is 0 up to 1996, then have this flat rate from 1996 onwards
  
  detec_array[1,,,27:31] <- linear_params[3]
  detec_array[2,,,27:31] <- linear_params[4]
  
  ## linear increase 2001 to 2009 ###
  
  detec_array[1,,,32:40] <- sweep(detec_array[1,,,32:40], 3, (0.02 + seq(1,9,1) * linear_params[5]), "+")
  detec_array[2,,,32:40] <- sweep(detec_array[2,,,32:40], 3, (0.02 + seq(1,9,1) * linear_params[6]), "+")
  detec_array[3:4,,,32:40] <- sweep(detec_array[3:4,,,32:40], 4, (detec_array[3,1,1,31] + seq(1,9,1) * linear_params[7]), "+")
  detec_array[5:7,,,32:40] <- sweep(detec_array[5:7,,,32:40], 4, (detec_array[5,1,1,31] + seq(1,9,1) * linear_params[8]), "+")
  
  ### linear increase 2010 to 2015 #####
  
  detec_array[1,,,41:46] <- sweep(detec_array[1,,,41:46], 3, (detec_array[1,1,1,40] + seq(1,6,1) * linear_params[9]), "+")
  detec_array[2,,,41:46] <- sweep(detec_array[2,,,41:46], 3, (detec_array[2,1,1,40] + seq(1,6,1) * linear_params[10]), "+")
  detec_array[3:4,,,41:46] <- sweep(detec_array[3:4,,,41:46], 4,(detec_array[3,1,1,40] + seq(1,6,1) * linear_params[11]), "+")
  detec_array[5:7,,,41:46] <- sweep(detec_array[5:7,,,41:46], 4,(detec_array[5,1,1,40] + seq(1,6,1) * linear_params[12]), "+")
  
  ### same values during proj period 2016 - 2021
  
  detec_array[1,,,47:52] <- detec_array[1,,,46] 
  detec_array[2,,,47:52] <- detec_array[2,,,46]
  detec_array[3:4,,,47:52] <- detec_array[3:4,,,46]
  detec_array[5:7,,,47:52] <- detec_array[5:7,,,46]
  
  fp$diagn_rate <- detec_array
  
  return(fp)
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
    nparam_incid <- 3
    fp$incidinput <- ilogistic(seq_len(fp$ss$PROJ_YEARS), exp(theta[1:2]), theta[3] - fp$ss$proj_start)
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

  if(fp$eppmod %in% c("logrw", "logrspline")){
    nparam_incid <- fp$numKnots + 1L
    beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = exp(as.vector(fp$rvec.spldes %*% beta)),
                  iota = transf_iota(theta[fp$numKnots+1], fp))
    fp[names(param)] <- param
  }

  if(fp$linear_diagnosis == "gamma"){
  fp$gamma_max <- exp(theta[nparam_incid+1])
  fp$delta_rate <- exp(theta[nparam_incid+2])

  fp <- cumgamma_diagn_rate(fp$gamma_max, fp$delta_rate, fp)
  }
  if(fp$linear_diagnosis == "12 param"){
    linear_params <- theta[nparam_incid+1:(length(theta)-nparam_incid)]
    
    fp <- cumulative_linear_diagn_rate(linear_params, fp) 
  }
  if(fp$linear_diagnosis == "knot_linear"){
    knot_params <- theta[nparam_incid+1:(length(theta) - nparam_incid)]
    
    fp <- knot_linear_diagn_rate(knot_params, fp)
  }
  
  fp
}

#' Log-liklihood for new diagnoses and AIDS deaths and ART inits if want to use that
ll_csavr <- function(theta, fp, likdat){
  
  fp <- create_param_csavr(theta, fp)

  mod <- simmod(fp)
  
  ll_aidsdeaths <- 0 
  
  if(fp$aidsdeath == T){

  mod_aidsdeaths <- colSums(attr(mod, "hivdeaths"),,2) ## getting some odd NAs in here sometimes, for the moment lets let NA = 0
  # mod_aidsdeaths[is.na(mod_aidsdeaths)] <- 0 
  
  ll_aidsdeaths <- with(likdat$aidsdeaths, sum(dpois(aidsdeaths, mod_aidsdeaths[idx] * (1 - prop_undercount), log=TRUE)))
  }
  
  ll_diagnoses <- 0
  
  if(fp$diagnoses_use == T){
  mod_diagnoses <- calc_diagnoses(mod, fp)
  
  if(fp$likelihood_cd4 == F){
  
  ll_diagnoses <- with(likdat$diagnoses, sum(dpois(diagnoses, mod_diagnoses[idx] * (1 - prop_undercount), log=TRUE)))
  
  }else{
    ## So what I'm doing here is first using cases when we have no cd4 counts, so I extract those years from the data and model
    ## Then we sum up the likelihood of these two vectors for the period before 2001 when no cd4 counts 
    
    ll_diagnoses <- 0
    pre_cd4_data <- likdat$diagnoses[likdat$diagnoses$year < fp$time_at_which_get_cd4_counts, ]
    pre_cd4_mod_data <- mod_diagnoses[mod_diagnoses$year < fp$time_at_which_get_cd4_counts, ]
    
    pre_cd4_mod_vals <- rep(0,(fp$time_at_which_get_cd4_counts - 1970))
    for(i in 1970:(fp$time_at_which_get_cd4_counts-1)){
      year_dat <- pre_cd4_mod_data[pre_cd4_mod_data$year == i,]
      year_diag <- sum(year_dat$diagnoses)
      pre_cd4_mod_vals[i-1969] <- year_diag
    }
    
    ## Also sometimes getting some weird NA/ NAN values in pre_cd4_mod_vals, will for the moment make these equal to 0
    # pre_cd4_mod_vals[is.nan(pre_cd4_mod_vals)] <- 0
    # pre_cd4_mod_vals[is.na(pre_cd4_mod_vals)] <- 0
    # 
    
    ll_diagnoses <- with(pre_cd4_data, sum(dpois(total_cases, pre_cd4_mod_vals[idx], log=TRUE)))

   ## Now from 2001 onwards we have cd4 counts, so below I extract the data from 2001 onwards from both data frames
   ## then I use the fact the model data only has 4 levels, and calculate the likelihood for each stage and add this 
   ## iteratively to the overall ll for our diagnoses data 
    
    post_cd4_mod_data <- mod_diagnoses[mod_diagnoses$year >= fp$time_at_which_get_cd4_counts, ]
    for(i in 1:fp$stages){
      current_stage <- as.character(i)
      likdat_stage <- likdat$diagnoses[likdat$diagnoses$year >= fp$time_at_which_get_cd4_counts, ]
      mod_stage <- post_cd4_mod_data[post_cd4_mod_data$stage == current_stage, ]
      index_to_check <- c(1:nrow(likdat_stage))
      
      col_index <- i + 3
      
      ll_diagnoses_current_stage <- sum(dpois(likdat_stage[,col_index], mod_stage$diagnoses[index_to_check], log = T))
      
      ll_diagnoses <- ll_diagnoses + ll_diagnoses_current_stage
    }
  }
  }
  ll_art_inits <- 0
  
  if(fp$artinit_use == T){
    
    mod_artinits <- calc_artinits(mod, fp)
  
    if(fp$likelihood_cd4 == F){
      ll_art_inits <- with(likdat$art_init, sum(dpois(total_art, mod_artinits[idx], log = T)))
    }else
  for(i in 1:fp$stages){
    current_stage <- as.character(i)
    mod_stage <- mod_artinits[mod_artinits$stage == current_stage, ]
    likdat_stage <- likdat$art_init[likdat$art_init$year >= 2006, c(1,2,i+2)]
    index_to_check <- likdat_stage$idx
    
    ll_art_inits_current_stage <- sum(dpois(likdat_stage[,3], mod_stage[index_to_check, 1], log = T))
    
    ll_art_inits <- ll_art_inits + ll_art_inits_current_stage
    
    }
  
  }
  
tot_likelihood <- ll_aidsdeaths + ll_diagnoses + ll_art_inits
  
return(tot_likelihood)

}

rw_mean <- rep(0.5, 52)
rw_sd <- rep(0.25, 52)

ilogistic_theta_mean <- c(-1, -10, 1995)
ilogistic_theta_sd <- c(5, 5, 10)

idbllogistic_theta_mean <- c(-1, -1, 1995, -10, -10)
idbllogistic_theta_sd <- c(5, 5, 10, 5, 5)

diagn_theta_mean <- c(3, -3)
diagn_theta_sd <- c(5, 5)

diagn_linear_theta_mean <- rep(0.5, 12)
diagn_linear_theta_sd <- rep(0.5,12)

diagn_knot_linear_mean <- c(0.1, 0.5, 3, 9, 15)
diagn_knot_linear_sd <- c(0.05, 0.2, 1, 2, 2.5)

logiota_pr_mean <- -13
logiota_pr_sd <- 5

sample_prior_csavr <- function(n, fp){

  mat_eppmod <- sample_prior_eppmod(n, fp)
  if(fp$linear_diagnosis == "12 param"){
    mat_diagn <- sample_prior_piecewise_diagn(n, fp)
  }
  if(fp$linear_diagnosis == "gamma"){
  mat_diagn <- sample_prior_diagn(n, fp)
  }
  if(fp$linear_diagnosis == "knot_linear"){
    mat_diagn <- sample_prior_knot_diagn(n, fp)
  }

  cbind(mat_eppmod, mat_diagn)
}

sample_prior_eppmod <- function(n, fp){
  
  if(fp$eppmod == "logrw"){
    nparam <- fp$numKnots + 1L

    mat <- matrix(NA, n, nparam)
    mat[,1] <- rnorm(n, 0.2, 1)  # u[1]
    mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, rw_prior_shape, rw_prior_rate)  # u[2:numKnots]
    mat[,fp$numKnots+1] <- sample_iota(n, fp)

  } else if(fp$eppmod == "logrspline") {

    nparam <- fp$numKnots + 1L

    mat <- matrix(NA, n, nparam)
    mat[,1] <- rnorm(n, 0.2, 1)  # u[1]
    mat[,2:fp$numKnots] <- bayes_rmvt(n, fp$numKnots-1,tau2_init_shape, tau2_init_rate)  # u[2:numKnots]
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
      theta_mean <- c(theta_mean, rw_mean)
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

sample_prior_piecewise_diagn <- function(n, fp){
  
  nparam <- length(diagn_linear_theta_mean)
  val <- rnorm(n * nparam, diagn_linear_theta_mean, diagn_linear_theta_sd)
  mat <- t(matrix(val, nparam, n))
  
  mat
}

sample_prior_knot_diagn <- function(n, fp){
  nparam <- length(diagn_knot_linear_mean)
  val <- rnorm(n * nparam, diagn_knot_linear_mean, diagn_knot_linear_sd)
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
  } else if(fp$eppmod == "logrspline"){
    lpr <- bayes_lmvt(theta_eppmod[rsp$fp$rtpenord:fp$numKnots], tau2_prior_shape, tau2_prior_rate)
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
  else if(fp$eppmod %in% c("logrw", "logrspline"))
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
