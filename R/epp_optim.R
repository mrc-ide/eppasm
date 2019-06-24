# extract from fitmod()
prepare_fp_for_fitmod <- function(epp, fp, likdat) {
    if(fp$eppmod %in% c("logrw", "rhybrid")) { # THIS IS REALLY MESSY, NEED TO REFACTOR CODE
      fp$SIM_YEARS <- as.integer(max(likdat$ancsite.dat$df$yidx,
                                     likdat$hhs.dat$yidx,
                                     likdat$ancrtcens.dat$yidx,
                                     likdat$hhsincid.dat$idx))
      fp$proj.steps <- seq(fp$ss$proj_start+0.5,
                           fp$ss$proj_start-1+fp$SIM_YEARS+0.5, 
                           by = 1/fp$ss$hiv_steps_per_year)
    } else {
      fp$SIM_YEARS <- fp$ss$PROJ_YEARS
    }

    ## Prepare the EPP model
    tsEpidemicStart <- ifelse(epp, fp$tsEpidemicStart, fp$ss$time_epi_start+0.5)
    if (fp$eppmod == "rspline")
      fp <- prepare_rspline_model(fp, tsEpidemicStart = tsEpidemicStart)
    else if (fp$eppmod == "rtrend")
      fp <- prepare_rtrend_model(fp)
    else if (fp$eppmod == "logrw")
      fp <- prepare_logrw(fp)
    else if (fp$eppmod == "rhybrid")
      fp <- prepare_rhybrid(fp)
    else if (fp$eppmod == "rlogistic")
      fp$tsEpidemicStart <-
        fp$proj.steps[which.min(abs(fp$proj.steps - fp$ss$time_epi_start+0.5))]

    fp$logitiota <- TRUE

    ## Prepare the incidence model
    fp$incidmod <- "eppspectrum"
    fp
}

# extract from fitmod()
epp_optim <- function(fp, likdat, opt_init, opt_method, opt_diffstep, opthess) {
  optfn <- function(theta, fp, likdat) {
    lprior(theta, fp) + sum(ll(theta, fp, likdat))
  }
  if (is.null(opt_init)) {
    X0       <- sample_prior(B0, fp)
    lpost0   <- likelihood(X0, fp, likdat, log=TRUE) + prior(X0, fp, log=TRUE)
    opt_init <- X0[which.max(lpost0)[1], ]
  }
  opt <- optim(opt_init, optfn, fp = fp, likdat = likdat, method = opt_method,
               control = list(fnscale = -1, trace = 4, maxit = opt_maxit,
                            ndeps = rep(opt_diffstep, length(opt_init))))
  opt$fp     <- fp
  opt$likdat <- likdat
  opt$param  <- fnCreateParam(opt$par, fp)
  opt$mod    <- simmod(update(fp, list=opt$param))
  if (opthess) {
    opt$hessian <- optimHess(opt_init, optfn, fp=fp, likdat=likdat,
                             control=list(fnscale=-1,
                                          trace=4,
                                          maxit=1e3,
                                          ndeps=rep(.Machine$double.eps^0.5,
                                                    length(opt_init))))
    opt$resample <- mvtnorm::rmvnorm(B.re, opt$par, solve(-opt$hessian))
  }
  return(opt)
}