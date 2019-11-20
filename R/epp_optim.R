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

obj_fn = function(theta, fp, likdat) 
  lprior(theta, fp) + sum(ll_all(theta, fp, likdat))

#' extract from fitmod()
#' 
#' fit model with optim()
#' 
#' @param epp ...
#' @param fp fix parameters
#' @param likdat data likelihood
#' @param control_optim a list specify as optim arguments, e.g. list(fn=,control=list())
#' @param B0 Number of parameter samples to find MAP as optim's starting value
#' @param B.re Number of resamples
#' @param doParallel use mclapply when finding initial value
epp_optim <- function(epp=FALSE, fp, likdat, control_optim, B0, B.re, doParallel) {
  .control.optim <- list(par=NULL, fn=obj_fn, method="BFGS", hessian=TRUE,
                         control=list(fnscale=-1, trace=4, maxit=1e3))
  if (!is.null(names(control_optim)))
    .control.optim <- modifyList(.control.optim, control_optim)
  if (is.null(.control.optim$par)) { # Find starting values that MAP
    X0     = eppasm:::sample.prior(B0, fp)
    eppasm:::likelihood(X0[1, ], fp, likdat, log=TRUE, doParallel)
    lpost0 = eppasm:::likelihood(X0, fp, likdat, log=TRUE, doParallel) + 
             eppasm:::prior(X0, fp, log=TRUE)
    .control.optim$par = X0[which.max(lpost0)[1], ]
    cat('best MAP', max(lpost0), '\n')
  }
  # .control.optim$ndeps <- rep(opt_diffstep, length(.control.optim$par)))
  .control.optim <- modifyList(.control.optim, list(fp = fp, likdat = likdat))
  opt = do.call("optim", .control.optim)
  opt$fp     = fp
  opt$likdat = likdat
  opt$param  = fnCreateParam(opt$par, fp)
  opt$mod    = simmod(update(fp, list=opt$param))
  optclass   = ifelse(epp, "eppopt", "specopt")
  if (.control.optim$hessian) {
    opt$resample = mvtnorm::rmvnorm(B.re, opt$par, solve(-opt$hessian))
    optclass     = c(optclass, ifelse(epp, "eppfit", "specfit"))
  }
  class(opt) = optclass
  opt
}