#' Prepare likelihood data inputs for CSAVR
#'
#' 
prepare_likdat_csavr <- function(csavrd, fp){
  
  csavrd$idx <- csavrd$year - fp$ss$proj_start + 1L

  # Patch to handle '_unreported' specification in <..MV5>
  if(!exists("aids_deaths_undercount", csavrd))
    csavrd$aids_deaths_undercount <- 0
  if(exists("aids_deaths_unreported", csavrd)){
    csavrd$aids_deaths_unreported[is.na(csavrd$aids_deaths_unreported)] <- 0
    csavrd$aids_deaths <- csavrd$aids_deaths + csavrd$aids_deaths_unreported
  }
  
  aidsdeaths <- setNames(data.frame(csavrd[c("year", "idx", "aids_deaths", "aids_deaths_undercount")]),
                         c("year", "idx", "aidsdeaths", "prop_undercount"))
  aidsdeaths <- subset(aidsdeaths, !is.na(aidsdeaths))
  aidsdeaths$prop_undercount[is.na(aidsdeaths$prop_undercount)] <- 0
  aidsdeaths$prop_undercount <- aidsdeaths$prop_undercount / 100

  if(fp$likelihood_cd4 == F){
  if(!exists("new_cases_undercount", csavrd))
    csavrd$new_cases_undercount <- 0
  if(exists("new_cases_unreported", csavrd)){
    csavrd$new_cases_unreported[is.na(csavrd$new_cases_unreported)] <- 0
    csavrd$new_cases <- csavrd$new_cases + csavrd$new_cases_unreported
  }

  diagnoses <- setNames(data.frame(csavrd[c("year", "idx", "total_cases", "new_cases_undercount")]),
                        c("year", "idx", "diagnoses", "prop_undercount"))
  diagnoses <- subset(diagnoses, !is.na(diagnoses))
  diagnoses$prop_undercount[is.na(diagnoses$prop_undercount)] <- 0
  diagnoses$prop_undercount <- diagnoses$prop_undercount / 100
  }else{
    diagnoses <- setNames(data.frame(csavrd[c("year","idx","total_cases",
                                              "stage_1_cases","stage_2_cases","stage_3_cases","stage_4_cases")]),
                          c("year","idx","total_cases","stage_1_cases","stage_2_cases","stage_3_cases","stage_4_cases"))
    
  }
  
  art_init <- NULL
  
  if(fp$artinit_use == T){
    ## So we now decide if we want to use all art_inits or cd4 art_inits 
    if(fp$likelihood ==F){
      art_init <- setNames(data.frame(csavrd[c("year","idx","total_art")]),
                           c("year","idx","total_art"))
      
    }else{
    
    art_init <- setNames(data.frame(csavrd[c("year","idx","stage_1_art","stage_2_art","stage_3_art","stage_4_art")]),
                         c("year","idx","stage_1_art","stage_2_art","stage_3_art","stage_4_art"))
    }
  }
  
  list(diagnoses = diagnoses, aidsdeaths = aidsdeaths, art_init = art_init)
}


#' Fit model
#'
#' @param ... Updates to fixed parameters (fp) object to specify fitting options.
fitmod_csavr <- function(obj, ..., B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500, opt_iter=0, 
                         sample_prior=eppasm:::sample.prior,
                         prior=eppasm:::prior,
                         likelihood=eppasm:::likelihood,
                         optfit=FALSE, opt_method="BFGS", opt_init=NULL, opt_maxit=1000, opt_diffstep=1e-3){

  ## Update fp
  fp <- update(obj$fp, ...)

  if(fp$eppmod == "logrw")
    fp <- prepare_logrw(fp, 1975.5)

  if(fp$eppmod == "logrspline")
    fp <- prepare_rspline_model(fp, tsEpidemicStart=1975.5)


  ## Prepare likelihood data
  likdat <- prepare_likdat_csavr(obj$csavrd, fp)

  ## Fit using optimization
  if(optfit){
    lpost <- function(theta, fp, likdat) lprior_csavr(theta, fp) + ll_csavr(theta, fp, likdat)
    if(is.null(opt_init)){
      X0 <- sample_prior_csavr(B0, fp)
      lpost0 <- likelihood_csavr(X0, fp, likdat, log=TRUE) + prior_csavr(X0, fp, log=TRUE)
      opt_init <- X0[which.max(lpost0)[1],]
    }
    start.time <- proc.time()
    fit <- optim(opt_init, lpost, fp=fp, likdat=likdat, method=opt_method, control=list(fnscale=-1, trace=4, maxit=opt_maxit, ndeps=rep(opt_diffstep, length(opt_init))))
    fit.time <- proc.time() - start.time
    fit$fp <- create_param_csavr(fit$par, fp)
    fit$mod <- simmod(fit$fp)
    class(fit) <- "csavropt"

  } else { # IMIS fit

    ## If IMIS fails, start again
    fit <- try(stop(""), TRUE)
    while(inherits(fit, "try-error")){
      start.time <- proc.time()
      fit <- try(imis(B0, B, B.re, number_k, opt_iter, fp=fp, likdat=likdat,
                      sample_prior = sample_prior_csavr,
                      prior = prior_csavr,
                      likelihood = likelihood_csavr,
                      dsamp = prior_csavr))
      fit.time <- proc.time() - start.time
    }
    fit$fp <- fp
    class(fit) <- "csavrfit"
  }

  fit$dat <- obj$csavrd
  fit$likdat <- likdat
  fit$time <- fit.time

  fit
}


calc_plhiv <- function(mod) colSums(mod[,,2,],,2)
calc_hivdeaths <- function(mod) colSums(attr(mod, "hivdeaths"),,2)
calc_rt <- function(mod, fp) incid(mod) / (prev(mod) * (1 - (1 - fp$relinfectART) * artcov15to49(mod)))

#' Simulate posterior model outputs
simfit.csavrfit <- function(fit){
  
  fp_list <- lapply(seq_len(nrow(fit$resample)), function(ii) create_param_csavr(fit$resample[ii,], fit$fp))
  mod_list <- lapply(fp_list, simmod)
  
  fit$rt <- mapply(calc_rt, mod = mod_list, fp = fp_list)
  fit$prev <- sapply(mod_list, prev)
  fit$incid <- mapply(incid, mod = mod_list, fp = fp_list)

  fit$diagnoses <- mapply(calc_diagnoses, mod = mod_list, fp = fp_list)
  idx <- fit$likdat$diagnoses$idx
  fit$reported_diagnoses <- array(NA, c(52,ncol(fit$diagnoses)))
  if(fp_list[[1]]$likelihood_cd4 == F){
  fit$reported_diagnoses[idx, ] <- fit$diagnoses[idx, ] * (1 - fit$likdat$diagnoses$prop_undercount)
  }else{
    for(i in 1:ncol(fit$diagnoses)){
      a <- fit$diagnoses[1, i] 
      b <- a$diagnoses
      tot_diag_across_cd4 <- b[1:52] + b[53:104] + b[105:156] + b[157:208]
      fit$reported_diagnoses[idx, i] <- tot_diag_across_cd4[idx] 
      
    }
  }
    
  fit$hivdeaths <- sapply(mod_list, calc_hivdeaths)
  idx <- fit$likdat$aidsdeaths$idx
  fit$reported_aidsdeaths <- array(NA, dim(fit$hivdeaths))
  fit$reported_aidsdeaths[idx, ] <- fit$hivdeaths[idx, ] * (1 - fit$likdat$aidsdeaths$prop_undercount)

  fit$plhiv <- sapply(mod_list, calc_plhiv)
  fit$popsize <- sapply(mod_list, colSums, dims=3)

  return(fit)
}

tidy.csavrfit <- function(fit){
  
  if(!exists("prev", fit))
    fit <- simfit(fit)
  
  years <- fit$fp$ss$proj_start + 1:fit$fp$ss$PROJ_YEARS - 1L

  vars <- c("outcome", "years", "mean", "se", "median", "lower", "upper",
            "lik_data", "vld_data")

  dfdiagn <- data.frame(outcome = "Reported diagnoses (15+)",
                        year = years,
                        estci2(fit$reported_diagnoses))
  if(fit$fp$likelihood_cd4 == F){
  dfdiagn <- merge(dfdiagn,
                   with(fit$likdat$diagnoses,
                        data.frame(year = year, lik_data = diagnoses)),
                   all.x=TRUE)
  }else{
    dfdiagn <- merge(dfdiagn,
                     with(fit$likdat$diagnoses,
                          data.frame(year = year, lik_data = total_cases)),
                     all.x=TRUE)
    
  }

  # dfplhiv <- data.frame(outcome = "PLHIV (15+)",
  #                       year = years,
  #                       estci2(fit$plhiv))
  # dfplhiv <- merge(dfplhiv,
  #                  setNames(fit$dat[c("year", "plhiv")],
  #                           c("year", "vld_data")))

  dfhivdeaths <- data.frame(outcome = "Reported AIDS deaths (15+)",
                            year = years,
                            estci2(fit$reported_aidsdeaths))
  dfhivdeaths <- merge(dfhivdeaths,
                       with(fit$likdat$aidsdeaths,
                            data.frame(year = year, lik_data = aidsdeaths)),
                       all.x=TRUE)

  dflogrt <- data.frame(outcome = "log r(t)",
                        year = years,
                        estci2(log(fit$rt)))
  
  dfincid <- data.frame(outcome = "Incidence (15-49)",
                        year = years,
                        estci2(fit$incid))
  
  dfprev <- data.frame(outcome = "Prevalence (15-49)",
                       year = years,
                       estci2(fit$prev))

                   
  df <- list(dfdiagn, dfhivdeaths, dflogrt, dfincid, dfprev)
  df <- lapply(df, function(x){x[setdiff(vars, names(x))] <- NA; x})
  
  do.call(rbind, df)
}

