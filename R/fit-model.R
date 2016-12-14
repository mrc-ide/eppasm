## Prepare national fit. Aggregates ANC data from regional EPP files.
prepare_national_fit <- function(pjnz, upd.path, proj.end=2013.5, hiv_steps_per_year = 10L){

  ## spectrum
  demp <- read_demog_param(upd.path)
  projp <- read_hivproj_param(pjnz)

  specfp <- create_spectrum_fixpar(projp, demp, proj_end = as.integer(proj.end), time_epi_start = projp$yr_start, hiv_steps_per_year= hiv_steps_per_year)  # Set time_epi_start tomatch EPP

  ## epp
  eppd <- read_epp_data(pjnz)
  epp.subp <- read_epp_subpops(pjnz)
  epp.input <- read_epp_input(pjnz)

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))
  val <- list()

  attr(val, "eppd") <- list(anc.used = do.call(c, lapply(eppd, "[[", "anc.used")),
                            anc.prev = do.call(rbind, lapply(eppd, "[[", "anc.prev")),
                            anc.n = do.call(rbind, lapply(eppd, "[[", "anc.n")))
  attr(val, "likdat") <- list(anclik.dat = with(attr(val, "eppd"), fnPrepareANCLikelihoodData(anc.prev, anc.n, anc.used, projp$yr_start)))
  attr(val, "likdat")$lastdata.idx <- max(unlist(attr(val, "likdat")$anclik.dat$anc.idx.lst),
                                          unlist(lapply(lapply(lapply(eppd, "[[", "hhs"), epp:::fnPrepareHHSLikData, projp$yr_start), "[[", "idx")))
  attr(val, "likdat")$firstdata.idx <- min(unlist(attr(val, "likdat")$anclik.dat$anc.idx.lst),
                                           unlist(lapply(lapply(lapply(eppd, "[[", "hhs"), epp:::fnPrepareHHSLikData, projp$yr_start), "[[", "idx")))
  attr(val, "specfp") <- specfp
  attr(val, "eppfp") <- fnCreateEPPFixPar(epp.input, proj.end = proj.end)
  attr(val, "country") <- attr(eppd, "country")

  return(val)
}


fitmod <- function(obj, ..., epp=FALSE, B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500, D=0, opt_iter=0,
                   sample.prior=eppspectrum:::sample.prior,
                   prior=eppspectrum:::prior,
                   likelihood=eppspectrum:::likelihood){

  ## ... : updates to fixed parameters (fp) object to specify fitting options

  if(epp)
    fp <- update(attr(obj, 'eppfp'), ...)
  else
    fp <- update(attr(obj, 'specfp'), ...)

  likdat <- attr(obj, 'likdat')


  ## If IMIS fails, start again
  fit <- try(stop(""), TRUE)
  while(inherits(fit, "try-error")){
    start.time <- proc.time()
    fit <- try(IMIS(B0, B, B.re, number_k, D, opt_iter, fp=fp, likdat=likdat,
                    sample.prior=sample.prior, prior=prior, likelihood=likelihood))
    fit.time <- proc.time() - start.time
  }
  fit$fp <- fp
  fit$likdat <- likdat
  fit$time <- fit.time

  class(fit) <- "specfit"

  return(fit)
}



## simulate incidence and prevalence
simfit.specfit <- function(fit, rwproj=FALSE, ageprevdat=TRUE, agegr3=TRUE, aidsdeaths=FALSE, pregprev=TRUE){
  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (fit$likdat$lastdata.idx-1)*fit$fp$ss$hiv_steps_per_year+1

    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- epp:::sim_rvec_rwproj(par$rvec, firstidx, lastidx, 1/fit$fp$ss$hiv_steps_per_year); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)
  
  fit$rvec <- sapply(mod.list, attr, "rvec")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp.list)
  fit$popsize <- sapply(mod.list, colSums, dims=3)

  if(pregprev)
    fit$pregprev <- sapply(mod.list, fnPregPrev)

  if(ageprevdat)
    fit$ageprevdat <- sapply(mod.list, ageprev, arridx=fit$likdat$hhsage.dat$arridx, agspan=5)

  if(agegr3){
    fit$agegr3prev <- lapply(mod.list, ageprev, aidx=c(15, 25, 35)-fit$fp$ss$AGE_START+1L, sidx=1:2,
                             yidx=(1999-fit$fp$ss$proj_start+1L):fit$fp$ss$PROJ_YEARS, agspan=c(10, 10, 15))
    fit$agegr3prev <- do.call(abind::abind, c(fit$agegr3prev, along=4))
  }

  if(aidsdeaths)
    fit$aidsdeaths <- sapply(lapply(mod.list, attr, "hivdeaths"), colSums, dims=2)
  
    
  return(fit)
}


sim_mod_list <- function(fit, rwproj=FALSE){

  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) eppspectrum::fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    dt <- if(inherits(fit$fp, "eppfp")) fit$fp$dt else 1.0/fit$fp$ss$hiv_steps_per_year
    
    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (fit$likdat$lastdata.idx-1)/dt+1

    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- epp:::sim_rvec_rwproj(par$rvec, firstidx, lastidx, dt); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)

  ## strip unneeded attributes to preserve memory
  mod.list <- lapply(mod.list, function(mod){ attributes(mod)[!names(attributes(mod)) %in% c("class", "dim", "infections", "hivdeaths", "natdeaths", "rvec")] <- NULL; mod})

  return(mod.list)
}

## ' aggregate lists of model fits 
aggr_specfit <- function(fitlist, rwproj=sapply(fitlist, function(x) x$fp$eppmod) == "rspline"){
  allmod <- parallel::mcmapply(sim_mod_list, fitlist, rwproj, SIMPLIFY=FALSE)

  modaggr <- lapply(do.call(mapply, c(FUN=list, allmod, SIMPLIFY=FALSE)), Reduce, f="+")
  ##
  infectionsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "infections"), SIMPLIFY=FALSE)), Reduce, f="+")
  hivdeathsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "hivdeaths"), SIMPLIFY=FALSE)), Reduce, f="+")
  natdeathsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "natdeaths"), SIMPLIFY=FALSE)), Reduce, f="+")
  ##
  modaggr <- mapply("attr<-", modaggr, "infections", infectionsaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "hivdeaths", hivdeathsaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "natdeaths", natdeathsaggr, SIMPLIFY=FALSE)
  ##
  modaggr <- mapply("attr<-", modaggr, "prev15to49", lapply(modaggr, calc_prev15to49, fitlist[[1]]$fp), SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "incid15to49", lapply(modaggr, calc_incid15to49, fitlist[[1]]$fp), SIMPLIFY=FALSE)
  ##
  modaggr <- lapply(modaggr, "class<-", c("specaggr", "spec"))
  return(modaggr)
}




####  Parameters  ####

frr_cd4_stage <- rbind(c(2.2092601, 1.0498458, 0.7909242, 0.7457322, 0.6370297, 0.5784465, 0.5784465, 0.5784465),
c(1.9957325, 0.9476719, 0.7137200, 0.6728340, 0.5747181, 0.5218241, 0.5218241, 0.5218241),
c(1.6978419, 0.8060810, 0.6068454, 0.5719130, 0.4884228, 0.4435285, 0.4435285, 0.4435285),
c(1.2546568, 0.5960523, 0.4484867, 0.4226725, 0.3609799, 0.3276454, 0.3276454, 0.3276454),
c(0.8321120, 0.3955791, 0.2974010, 0.2802773, 0.2391743, 0.2170073, 0.2170073, 0.2170073),
c(0.6498167, 0.3092457, 0.2323928, 0.2188925, 0.1869424, 0.1695819, 0.1695819, 0.1695819),
c(0.6196790, 0.2947683, 0.2214894, 0.2086464, 0.1781029, 0.1616417, 0.1616417, 0.1616417))

frr_art_stage <- array(0, c(3,7,8))
frr_art_stage[1,,] <- frr_cd4_stage
frr_art_stage[2,,] <- frr_cd4_stage
frr_art_stage[3,,] <- 0.8


frr_cd4_stage1 <- frr_cd4_stage
frr_art_stage1 <- frr_art_stage
frr_art_stage1[3,,] <- 1.0

frr_cd4_default <- rbind(c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47),
                         c(1.2, 1.2, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47))

frr_art_default <- array(0, c(3,7,8))
frr_art_default[1,,] <- frr_cd4_default
frr_art_default[2,,] <- frr_cd4_default
frr_art_default[3,,] <- 1.0
