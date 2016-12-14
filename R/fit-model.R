## Prepare fit by EPP regions
prepare_spec_fit <- function(pjnz, proj.end=2016.5){

  ## epp
  eppd <- read_epp_data(pjnz)
  epp.subp <- read_epp_subpops(pjnz)
  epp.input <- read_epp_input(pjnz)

  epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)

  ## spectrum
  demp <- read_specdp_demog_param(pjnz)
  projp <- read_hivproj_param(pjnz)

  specfp.subp <- create_subpop_specfp(projp, demp, eppd, proj_end=proj.end)
  

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat, anchor.year=epp.input$start.year))
  val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "specfp", specfp.subp)
  val <- set.list.attr(val, "country", attr(eppd, "country"))
  val <- set.list.attr(val, "region", names(eppd))

  return(val)
}


create_subpop_specfp <- function(projp, demp, eppd, ..., popadjust=TRUE){

  ## Urban/Rural target population database
  urpop <- readRDS(system.file("extdata", "urpop.RDS", package="eppspectrum"))
  
  country <- attr(eppd, "country")

  ## Update demp for subpopulation 
  demp.subpop <- list()
  for(subpop in names(eppd)){
    subp <- if(subpop %in% c("Urbain", "Urbaine", "Urban")) "Urban"
            else if(subpop %in%  c("Rural", "Rurale")) "Rural"
            else subpop  # bloody French...
    demp.subpop[[subpop]] <- demp
    demp.subpop[[subpop]]$basepop <- urpop[,,,subp, country]
    demp.subpop[[subpop]]$netmigr[] <- 0
  }

  ## Apportion ART
  ## If national survey data are available, apportion ART according to relative average HH survey prevalence in each subpopulation,
  ## If no HH survey, apportion based on relative mean ANC prevalence

  get15to49pop <- function(demp, year) sum(demp$basepop[as.character(15:49),,as.character(year)])
  subpop.dist <- prop.table(sapply(demp.subpop, get15to49pop, 2010))
  
  if(nrow(subset(eppd[[1]]$hhs, used)) != 0){ # HH survey data available
    hhsprev.means <- sapply(lapply(eppd, function(dat) na.omit(dat$hhs$prev[dat$hhs$used])), mean)
    art.dist <- prop.table(subpop.dist * hhsprev.means)
  } else {  ## no HH survey data
    ## Apportion ART according to relative average ANC prevalence in each subpopulation
    ancprev.means <- sapply(lapply(eppd, "[[", "anc.prev"), mean, na.rm=TRUE)
    art.dist <- prop.table(subpop.dist * ancprev.means)
  }

  ## Update projp for subpopulation
  projp.subpop <- list()
  for(subpop in names(eppd)){
    projp.subpop[[subpop]] <- projp
    isartnum <- projp$art15plus_numperc == 0
    projp.subpop[[subpop]]$art15plus_num[isartnum] <- projp$art15plus_num[isartnum] * art.dist[subpop]
  }

  specfp.subpop <- list()
  for(subpop in names(eppd))
    specfp.subpop[[subpop]] <- create_spectrum_fixpar(projp.subpop[[subpop]], demp.subpop[[subpop]], ..., popadjust=popadjust)

  return(specfp.subpop)
}


## Prepare national fit. Aggregates ANC data from regional EPP files.
prepare_national_fit <- function(pjnz, upd.path=NULL, proj.end=2013.5, hiv_steps_per_year = 10L){

  ## spectrum
  if(!is.null(upd.path))
    demp <- read_demog_param(upd.path)
  else
    demp <- read_specdp_demog_param(pjnz)
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
  attr(val, "likdat") <- list(anclik.dat = with(attr(val, "eppd"), anclik::fnPrepareANCLikelihoodData(anc.prev, anc.n, anc.used, projp$yr_start)))
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
                   sample.prior=epp:::sample.prior,
                   prior=epp:::prior,
                   likelihood=epp:::likelihood){

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

  if(epp)
    class(fit) <- "eppfit"
  else
    class(fit) <- "specfit"

  return(fit)
}



## simulate incidence and prevalence
simfit.specfit <- function(fit, rwproj=fit$fp$eppmod == "rspline", ageprevdat=FALSE, agegr3=FALSE, aidsdeaths=FALSE, pregprev=TRUE){
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


sim_mod_list <- function(fit, rwproj=fit$fp$eppmod == "rspline"){

  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

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
  mod.list <- lapply(mod.list, function(mod){ attributes(mod)[!names(attributes(mod)) %in% c("class", "dim", "infections", "hivdeaths", "natdeaths", "rvec", "popadjust")] <- NULL; mod})

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

calc_prev15to49 <- function(mod, fp){
  colSums(mod[fp$ss$p.age15to49.idx,,2,],,2)/colSums(mod[fp$ss$p.age15to49.idx,,,],,3)
}

calc_incid15to49 <- function(mod, fp){
  c(0, colSums(attr(mod, "infections")[fp$ss$p.age15to49.idx,,-1],,2)/colSums(mod[fp$ss$p.age15to49.idx,,1,-fp$ss$PROJ_YEARS],,2))
}

calc_pregprev <- function(mod, fp){
  warning("not yet implemented")
}
