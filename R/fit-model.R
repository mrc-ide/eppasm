#' Prepare fit by EPP regions
#'
#' @param pjnz file path to Spectrum PJNZ file.
#' @param proj.end end year for projection.
#' @param popadjust popadjust data
#' @param popupdate logical should target population be updated to match age-specific population size from DP file and Urban from EPP XML.
#' @param use_ep5 see `read_specdp_demog_param`
prepare_spec_fit <- function(pjnz, proj.end=2016.5, popadjust = NULL, popupdate=TRUE, use_ep5=FALSE){

  ## epp
  eppd <- epp::read_epp_data(pjnz)
  ## epp.subp <- epp::read_epp_subpops(pjnz)
  ## epp.input <- epp::read_epp_input(pjnz)

  ## epp.subp.input <- epp::fnCreateEPPSubpops(epp.input, epp.subp, eppd)

  country <- attr(eppd, "country")
  cc <- attr(eppd, "country_code")

  ## melt site-level data
  eppd <- Map("[[<-", eppd, "ancsitedat", lapply(eppd, melt_ancsite_data))

  ## tidy HHS data
  eppd <- Map("[[<-", eppd, "hhs", lapply(eppd, tidy_hhs_data))

  attr(eppd, "country") <- country
  attr(eppd, "country_code") <- cc
    
  ## spectrum
  demp <- read_specdp_demog_param(pjnz, use_ep5=use_ep5)
  projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)
  epp_t0 <- read_epp_t0(pjnz)

  ## If popadjust = NULL, look for subp if more than 1 EPP region
  if(is.null(popadjust))
    popadjust <- length(eppd) > 1

  ## If Urban/Rural fit, read percentage urban from EPP XML file
  if(length(eppd) == 2 && all(sort(substr(names(eppd), 1, 1)) == c("R", "U")))
    perc_urban <- read_epp_perc_urban(pjnz) / 100
  else
    perc_urban <- NULL
    
  specfp.subp <- create_subpop_specfp(projp, demp, eppd, proj_end=proj.end, epp_t0=epp_t0,
                                      popadjust = popadjust, popupdate = popupdate, perc_urban = perc_urban)
  

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  ## val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, epp::fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "specfp", specfp.subp)
  val <- set.list.attr(val, "country", read_country(pjnz))
  val <- set.list.attr(val, "region", names(eppd))

  attr(val, "country") <- read_country(pjnz)
  attr(val, "region") <- read_region(pjnz)

  return(val)
}

#' Melt ANC-SS and site-level ANC-RT to long dataset
#' 
#' @param eppd EPP data
melt_ancsite_data <- function(eppd){
  
  anc.used <- data.frame(site=rownames(eppd$anc.prev), used=eppd$anc.used)
  anc.prev <- subset(reshape2::melt(eppd$anc.prev, varnames=c("site", "year"), value.name="prev"), !is.na(prev))
  anc.n <- subset(reshape2::melt(eppd$anc.n, varnames=c("site", "year"), value.name="n"), !is.na(n))
  
  ancsitedat <- merge(anc.used, anc.prev)
  ancsitedat <- merge(ancsitedat, anc.n)
  ancsitedat$type <- "ancss"

  if(exists("ancrtsite.prev", eppd) && !is.null(eppd$ancrtsite.prev)){
    ancrtsite.prev <- subset(reshape2::melt(eppd$ancrtsite.prev, varnames=c("site", "year"), value.name="prev"), !is.na(prev))
    ancrtsite.n <- subset(reshape2::melt(eppd$ancrtsite.n, varnames=c("site", "year"), value.name="n"), !is.na(n))
  
    ancrtsite <- merge(anc.used, ancrtsite.prev)
    ancrtsite <- merge(ancrtsite, ancrtsite.n)
    ancrtsite$type <- rep("ancrt", nrow(ancrtsite))
    
    ancsitedat <- rbind(ancsitedat, ancrtsite)
  }

  ancsitedat <- subset(ancsitedat, used)
  ancsitedat$agegr <- rep("15-49", nrow(ancsitedat))
  ancsitedat$age <- rep(15, nrow(ancsitedat))
  ancsitedat$agspan <- rep(35, nrow(ancsitedat))

  ancsitedat
}

tidy_hhs_data <- function(eppd){

  hhs <- eppd$hhs

  hhs$deff <- hhs$deff_approx <- rep(2.0, nrow(hhs)) # irrelevant assumption 
  hhs$n <- hhs$deff * hhs$prev * (1-hhs$prev) / hhs$se^2
  hhs$agegr <- rep("15-49", nrow(hhs))
  hhs$sex <- rep("both", nrow(hhs))
  
  hhs <- hhs[c("year", "sex", "agegr", "n", "prev", "se", "deff", "deff_approx", "used")]
  
  hhs
}

create_subpop_specfp <- function(projp, demp, eppd, epp_t0=setNames(rep(1975, length(eppd)), names(eppd)), ..., popadjust=TRUE, popupdate=TRUE, perc_urban=NULL){

  country <- attr(eppd, "country")
  country_code <- attr(eppd, "country_code")

  ## Update demp for subpopulation 
  demp.subpop <- list()
  gfr_subpop <- list()
  for(subpop in names(eppd)){
    ## if(country != "Malawi")
    strsubp <- if(subpop %in% c("Urbain", "Urbaine", "Urban")) "U"
               else if(subpop %in%  c("Rural", "Rurale")) "R"
               else subpop  # bloody French...
    demp.subpop[[subpop]] <- demp
    if (popadjust) {
      demp.subpop[[subpop]]$basepop <- subp[[grep(paste0("^", country_code, "_"), names(subp))]][[strsubp]][,,dimnames(demp$basepop)[[3]]]
      demp.subpop[[subpop]]$netmigr[] <- 0

      ## Record GFR for each subpop
      gfr_subpop[[subpop]] <- subset(subp_gfr, cc == country_code & eppregion == strsubp, c(gfr, survyear))
    }
  }

  ## Compare demp population with subpopulations
  ann_aggr_subpop <- colSums(Reduce("+", lapply(demp.subpop, "[[", "basepop")),,2)
  ann_demp_bp <- colSums(demp$basepop,,2)
  if(any(abs(ann_aggr_subpop[names(ann_demp_bp)] / ann_demp_bp - 1.0) > 0.05))
    warning(paste("National popultion population differs from aggregated subpopulations >5%:", country))


  if(popupdate){
    ## Rake subpopulation size to match DP population by age and sex, and area
    ## population distribution.

    if(!all(sapply(lapply(lapply(demp.subpop, "[[", "basepop"), dim),
                   identical, dim(demp$basepop))))
      stop(paste("Dimensions of basepop do not match subpopulation dimensions:",
                 attr(eppd, "country")))

    subpops <- do.call(abind::abind, c(lapply(demp.subpop, "[[", "basepop"), along=4))
    
    ## adjusted population sizes
    if(!is.null(perc_urban))
      areapop <- colSums(demp$basepop,,2) * cbind(Urban=perc_urban, Rural=1-perc_urban)
    else
      areapop <- colSums(demp$basepop,,2) * prop.table(colSums(subpops,,2), 1)        
    agesexpop <- demp$basepop

    ## Iteratively rescale population until difference < 0.1%
    while(any(abs(rowSums(subpops,,3) / agesexpop - 1.0) > 0.001, na.rm=TRUE)){
      
      ## Scale supopulation size to match national population by age/sex
      subpops <- subpops <- sweep(subpops, 1:3, agesexpop / rowSums(subpops,,3), "*")
      ## Scale subpopulations to match subpopulation distribution
      subpops <- sweep(subpops, 3:4, areapop / colSums(subpops,,2), "*")
    }

    for(subpop in names(demp.subpop))
      demp.subpop[[subpop]]$basepop <- subpops[,,,subpop]
  }

  if(length(demp.subpop) > 1){
    ## Apportion births according to population size and GFR
    for(subpop in names(demp.subpop)){
      survyear <- gfr_subpop[[subpop]]$survyear
      gfr <- gfr_subpop[[subpop]]$gfr
      fpop15to44 <- sum(demp.subpop[[subpop]]$basepop[as.character(15:44), "Female", as.character(survyear)])
      prop_births <- fpop15to44*gfr / demp$births[as.character(survyear)]
      
      demp.subpop[[subpop]]$births <- prop_births * demp$births
    }
    
    ## Rake births to national births
    births_adjust <- demp$births / rowSums(sapply(demp.subpop, "[[", "births"))
    for(subpop in names(demp.subpop))
      demp.subpop[[subpop]]$births <- births_adjust * demp.subpop[[subpop]]$births 
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

  ## Apportion age 14 HIV population
  ## Allocate relative to HIV prevalence and population size, same as ART population
  for(subpop in names(eppd)){
    projp.subpop[[subpop]]$age15hivpop <- projp.subpop[[subpop]]$age15hivpop * art.dist[subpop]
  }
  
  specfp.subpop <- list()
  for(subpop in names(eppd))
    specfp.subpop[[subpop]] <- create_spectrum_fixpar(projp.subpop[[subpop]], demp.subpop[[subpop]], ..., popadjust=popadjust, time_epi_start=epp_t0[subpop])

  return(specfp.subpop)
}


## Prepare national fit. Aggregates ANC data from regional EPP files.
prepare_national_fit <- function(pjnz, upd.path=NULL, proj.end=2013.5, hiv_steps_per_year = 10L, use_ep5=FALSE){

  ## spectrum
  if(!is.null(upd.path))
    demp <- read_demog_param(upd.path)
  else
    demp <- read_specdp_demog_param(pjnz, use_ep5=use_ep5)
  projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)
  epp_t0 <- read_epp_t0(pjnz)

  specfp <- create_spectrum_fixpar(projp, demp, proj_end = as.integer(proj.end), time_epi_start = epp_t0[1], hiv_steps_per_year= hiv_steps_per_year)  # Set time_epi_start to match first EPP population

  ## epp
  eppd <- epp::read_epp_data(pjnz)
  epp.subp <- epp::read_epp_subpops(pjnz)
  epp.input <- epp::read_epp_input(pjnz)

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))
  val <- list()

  ## aggregate census data across regions
  ancrtcens <- do.call(rbind, lapply(eppd, "[[", "ancrtcens"))
  if(!is.null(ancrtcens)) {
    ancrtcens <- subset(ancrtcens, !is.na(prev) & !is.na(n))
    if(nrow(ancrtcens)){
      ancrtcens$x <- ancrtcens$prev * ancrtcens$n
      ancrtcens <- aggregate(cbind(x,n) ~ year, ancrtcens, sum)
      ancrtcens$prev <- ancrtcens$x / ancrtcens$n
      ancrtcens <- ancrtcens[c("year", "prev", "n")]
    }
  }
  
  attr(val, "eppd") <- list(anc.used = do.call(c, lapply(eppd, "[[", "anc.used")),
                            anc.prev = do.call(rbind, lapply(eppd, "[[", "anc.prev")),
                            anc.n = do.call(rbind, lapply(eppd, "[[", "anc.n")),
                            ancrtsite.prev = do.call(rbind, lapply(eppd, "[[", "ancrtsite.prev")),
                            ancrtsite.n = do.call(rbind, lapply(eppd, "[[", "ancrtsite.n")),
                            ancrtcens = ancrtcens)

  attr(val, "specfp") <- specfp
  attr(val, "eppfp") <- epp::fnCreateEPPFixPar(epp.input, proj.end = proj.end)
  attr(val, "country") <- attr(eppd, "country")

  return(val)
}


fitmod <- function(obj, ..., epp=FALSE, B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500, opt_iter=0, 
                   sample_prior=eppasm:::sample.prior,
                   prior=eppasm:::prior,
                   likelihood=eppasm:::likelihood,
                   optfit=FALSE, opt_method="BFGS", opt_init=NULL, opt_maxit=1000, opt_diffstep=1e-3, opthess=TRUE){

  ## ... : updates to fixed parameters (fp) object to specify fitting options

  if(epp)
    fp <- update(attr(obj, 'eppfp'), ...)
  else
    fp <- update(attr(obj, 'specfp'), ...)


  ## Prepare likelihood data
  eppd <- attr(obj, "eppd")

  has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
  has_ancrtcens <- !is.null(eppd$ancrtcens) && nrow(eppd$ancrtcens)
  
  if(!has_ancrtsite)
    fp$ancrtsite.beta <- 0

  if(has_ancrtsite & has_ancrtcens)
    fp$ancrt <- "both"
  else if(has_ancrtsite & !has_ancrtcens)
    fp$ancrt <- "site"
  else if(!has_ancrtsite & has_ancrtcens)
    fp$ancrt <- "census"
  else
    fp$ancrt <- "none"

  if(epp)
    eppd$hhsage <- eppd$sibmx <- NULL

  likdat <- prepare_likdat(eppd, fp)
  fp$ancsitedata <- as.logical(nrow(likdat$ancsite.dat$df))

  if(fp$eppmod %in% c("logrw", "rhybrid")) { # THIS IS REALLY MESSY, NEED TO REFACTOR CODE
    
    fp$SIM_YEARS <- as.integer(max(likdat$ancsite.dat$df$yidx,
                                   likdat$hhs.dat$yidx,
                                   likdat$ancrtcens.dat$yidx,
                                   likdat$hhsincid.dat$idx,
                                   likdat$sibmx.dat$idx))
      
    fp$proj.steps <- seq(fp$ss$proj_start+0.5, fp$ss$proj_start-1+fp$SIM_YEARS+0.5, by=1/fp$ss$hiv_steps_per_year)
  } else
    fp$SIM_YEARS <- fp$ss$PROJ_YEARS

  ## Prepare the EPP model
  tsEpidemicStart <- if(epp) fp$tsEpidemicStart else fp$ss$time_epi_start+0.5
  if(fp$eppmod == "rspline")
    fp <- prepare_rspline_model(fp, tsEpidemicStart=tsEpidemicStart)
  else if(fp$eppmod == "rtrend")
    fp <- prepare_rtrend_model(fp)
  else if(fp$eppmod == "logrw")
    fp <- prepare_logrw(fp)
  else if(fp$eppmod == "rhybrid")
    fp <- prepare_rhybrid(fp)
  else if(fp$eppmod == "rlogistic")
    fp$tsEpidemicStart <- fp$proj.steps[which.min(abs(fp$proj.steps - fp$ss$time_epi_start+0.5))]

  fp$logitiota <- TRUE

  ## Prepare the incidence model
  if(exists("incidmod", where=fp) && fp$incidmod == "transm"){
    if(!exists("relsexact_cd4cat", where=fp))
      fp$relsexact_cd4cat <- c(1.0, 0.92, 0.76, 0.76, 0.55, 0.55, 0.55)
  } else
    fp$incidmod <- "eppspectrum"

  fp <- prepare_irr_model(fp)

  ## Fit using optimization
  if(optfit){
    optfn <- function(theta, fp, likdat) lprior(theta, fp) + sum(ll(theta, fp, likdat))
    if(is.null(opt_init)){
      X0 <- sample_prior(B0, fp)
      lpost0 <- likelihood(X0, fp, likdat, log=TRUE) + prior(X0, fp, log=TRUE)
      opt_init <- X0[which.max(lpost0)[1],]
    }
    opt <- optim(opt_init, optfn, fp=fp, likdat=likdat, method=opt_method, control=list(fnscale=-1, trace=4, maxit=opt_maxit, ndeps=rep(opt_diffstep, length(opt_init))))
    opt$fp <- fp
    opt$likdat <- likdat
    opt$param <- fnCreateParam(opt$par, fp)
    opt$mod <- simmod(update(fp, list=opt$param))
    if(opthess){
      opt$hessian <- optimHess(opt_init, optfn, fp=fp, likdat=likdat,
                               control=list(fnscale=-1,
                                            trace=4,
                                            maxit=1e3,
                                            ndeps=rep(.Machine$double.eps^0.5, length(opt_init))))
      opt$resample <- mvtnorm::rmvnorm(B.re, opt$par, solve(-opt$hessian))
    }

    if(epp)
      class(opt) <- "eppopt"
    else
      class(opt) <- "specopt"

    if(opthess){
      if(epp)
        class(opt) <- c(class, "eppfit")
      else
        class(opt) <- c(class, "specfit")
    }

    return(opt)
  }

  ## If IMIS fails, start again
  fit <- try(stop(""), TRUE)
  while(inherits(fit, "try-error")){
    start.time <- proc.time()
    fit <- try(imis(B0, B, B.re, number_k, opt_iter, fp=fp, likdat=likdat,
                    sample_prior=sample.prior, prior=prior, likelihood=likelihood))
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

get_likdat_range <- function(likdat) {

  firstdata.idx <- as.integer(min(likdat$ancsite.dat$df$yidx,
                                  likdat$hhs.dat$yidx,
                                  likdat$ancrtcens.dat$yidx,
                                  likdat$hhsincid.dat$idx,
                                  likdat$sibmx.dat$idx))

  lastdata.idx <- as.integer(max(likdat$ancsite.dat$df$yidx,
                                 likdat$hhs.dat$yidx,
                                 likdat$ancrtcens.dat$yidx,
                                 likdat$hhsincid.dat$idx,
                                 likdat$sibmx.dat$idx))

  c(firstdata.idx, lastdata.idx)
}

rw_projection <- function(fit) {

  if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
    stop("Random-walk projection is only used with r-spline model")
  
  datrange <- get_likdat_range(fit$likdat)
  firstidx <- datrange[1] * fit$fp$ss$hiv_steps_per_year
      lastidx <- (datrange[2]-1)*fit$fp$ss$hiv_steps_per_year+1
  
  ## replace rvec with random-walk simulated rvec
  fit$param <- lapply(fit$param, function(par){par$rvec <- epp:::sim_rvec_rwproj(par$rvec, firstidx, lastidx, 1/fit$fp$ss$hiv_steps_per_year); par})

  fit
}


## simulate incidence and prevalence
#' @importFrom abind abind
simfit.specfit <- function(fit,
                           rwproj=fit$fp$eppmod == "rspline",
                           ageprevdat=FALSE,
                           agegr3=FALSE,
                           artcov=FALSE,
                           ancartcov=FALSE,
                           mxoutputs=FALSE,
                           aidsdeaths=FALSE,
                           pregprev=TRUE,
                           ageincid=TRUE,
                           ageinfections=TRUE,
                           relincid = TRUE,
                           entrantprev=TRUE,
                           mod.list=NULL){

  if(is.null(mod.list)){
    fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

    
    if(rwproj)
      fit <- rw_projection(fit)

    fp_list <- lapply(fit$param, function(par) update(fit$fp, list=par))
    mod.list <- lapply(fp_list, simmod)

  } else {
    fp_list <- rep(fit$fp, length(mod.list))
  }
  
  fit$rvec <- sapply(mod.list, attr, "rvec_ts")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp_list)
  fit$popsize <- sapply(mod.list, colSums, dims=3)

  if(inherits(pregprev, "data.frame"))
    fit$pregprev <- mapply(agepregprev, mod=mod.list, fp=fp_list,
                           MoreArgs=list(aidx=pregprev$aidx, yidx=pregprev$yidx, agspan=pregprev$agspan))
  else if(is.logical(pregprev) && pregprev == TRUE)

  if(entrantprev)
    fit$entrantprev <- sapply(mod.list, attr, "entrantprev")

  if(is.logical(ageprevdat) && ageprevdat == TRUE)
    ageprevdat <- fit$likdat$hhs.dat
  if(inherits(ageprevdat, "data.frame")){
    fit$ageprevdat <- mapply(ageprev, mod=mod.list,
                             MoreArgs=list(aidx=ageprevdat$aidx,
                                           sidx=ageprevdat$sidx,
                                           yidx=ageprevdat$yidx,
                                           agspan=ageprevdat$agspan))
    fit$ageprevdat <- matrix(fit$ageprevdat, nrow=nrow(ageprevdat))
  }
  
  if(inherits(agegr3, "data.frame"))
    fit$agegr3prev <- mapply(ageprev, mod=mod.list,
                             MoreArgs=list(aidx=agegr3$aidx,
                                           sidx=agegr3$sidx,
                                           yidx=agegr3$yidx,
                                           agspan=agegr3$agspan))
  else if(is.logical(agegr3) && agegr3 == TRUE) {
    fit$agegr3prev <- lapply(mod.list, ageprev, aidx=c(15, 25, 35)-fit$fp$ss$AGE_START+1L, sidx=1:2,
                             yidx=(1999-fit$fp$ss$proj_start+1L):fit$fp$ss$PROJ_YEARS, agspan=c(10, 10, 15), expand=TRUE)
    fit$agegr3prev <- do.call(abind::abind, c(fit$agegr3prev, along=4))
  }

  if(aidsdeaths)
    fit$aidsdeaths <- sapply(lapply(mod.list, attr, "hivdeaths"), colSums, dims=2)

  if(artcov)
    fit$artcov <- sapply(mod.list, artcov15plus)

  if(ancartcov)
    fit$ancartcov <- mapply(agepregartcov, mod.list, fp_list, MoreArgs=list(aidx=1, yidx=1:fit$fp$ss$PROJ_YEARS, agspan=35, expand=TRUE))

  if(mxoutputs){
    fit$agemx <- abind::abind(lapply(mod.list, agemx), rev.along=0)
    dimnames(fit$agemx) <- with(fit$fp$ss, list(age=with(fit$fp$ss, AGE_START+0:(pAG-1)), sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
    
    fit$natagemx <- abind::abind(lapply(mod.list, agemx, nonhiv=TRUE), rev.along=0)
    dimnames(fit$natagemx) <- with(fit$fp$ss, list(age=with(fit$fp$ss, AGE_START+0:(pAG-1)), sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
    
    fit$q4515 <- abind::abind(lapply(mod.list, calc_nqx, fp=fit$fp, n=45, x=15), rev.along=0)
    dimnames(fit$q4515) <- with(fit$fp$ss, list(sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
    
    fit$q3515 <- abind::abind(lapply(mod.list, calc_nqx, fp=fit$fp, n=35, x=15), rev.along=0)
    dimnames(fit$q3515) <- with(fit$fp$ss, list(sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
    
    fit$natq4515 <- abind::abind(lapply(mod.list, calc_nqx, fp=fit$fp, n=45, x=15, nonhiv=TRUE), rev.along=0)
    dimnames(fit$natq4515) <- with(fit$fp$ss, list(sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
    
    fit$natq3515 <- abind::abind(lapply(mod.list, calc_nqx, fp=fit$fp, n=35, x=15, nonhiv=TRUE), rev.along=0)
    dimnames(fit$natq3515) <- with(fit$fp$ss, list(sex=c("male", "female"), proj_start + 0:(PROJ_YEARS-1), NULL))
  }


  if(ageincid){

    startyr <- fit$fp$ss$proj_start
    endyr <- fit$fp$ss$proj_start + fit$fp$ss$PROJ_YEARS-1L

    agegr3 <- lapply(mod.list, ageincid, aidx=c(15, 25, 35, 50)-fit$fp$ss$AGE_START+1L, sidx=1:2,
                     yidx=1999:endyr - startyr+1L, agspan=c(10, 10, 15, 31))
    age15to49 <- lapply(mod.list, ageincid, aidx=15-fit$fp$ss$AGE_START+1L, sidx=1:2,
                        yidx=1999:endyr - startyr+1L, agspan=35)
    age15plus <- lapply(mod.list, ageincid, aidx=15-fit$fp$ss$AGE_START+1L, sidx=1:2,
                        yidx=1999:endyr - startyr+1L, agspan=66)
    agegr3 <- abind::abind(agegr3, rev.along=0)
    age15to49 <- abind::abind(age15to49, rev.along=0)
    age15plus <- abind::abind(age15plus, rev.along=0)
    fit$agegr3incid <- abind::abind(agegr3, age15to49, age15plus, along=1)

    dimnames(fit$agegr3incid)[1:3] <- list(agegr=c("15-24", "25-34", "35-49", "50+", "15-49", "15+"),
                                           sex=c("male", "female"),
                                           year=1999:endyr)
  }

  if(ageinfections){
    
    startyr <- fit$fp$ss$proj_start
    endyr <- fit$fp$ss$proj_start + fit$fp$ss$PROJ_YEARS-1L
    
    agegr3 <- lapply(mod.list, ageinfections, aidx=c(15, 25, 35, 50)-fit$fp$ss$AGE_START+1L, sidx=1:2,
                     yidx=1999:endyr - startyr+1L, agspan=c(10, 10, 15, 31))
    age15to49 <- lapply(mod.list, ageinfections, aidx=15-fit$fp$ss$AGE_START+1L, sidx=1:2,
                        yidx=1999:endyr - startyr+1L, agspan=35)
    age15plus <- lapply(mod.list, ageinfections, aidx=15-fit$fp$ss$AGE_START+1L, sidx=1:2,
                        yidx=1999:endyr - startyr+1L, agspan=66)
    agegr3 <- abind::abind(agegr3, rev.along=0)
    age15to49 <- abind::abind(age15to49, rev.along=0)
    age15plus <- abind::abind(age15plus, rev.along=0)
    fit$agegr3infections <- abind::abind(agegr3, age15to49, age15plus, along=1)

    dimnames(fit$agegr3infections)[1:3] <- list(agegr=c("15-24", "25-34", "35-49", "50+", "15-49", "15+"),
                                           sex=c("male", "female"),
                                           year=1999:endyr)
  }

  if(relincid){

    ## Incidence relative to 25-29y
    fit$relincid <- lapply(mod.list, ageincid, aidx=3:9*5-fit$fp$ss$AGE_START+1L, sidx=1:2,
                           yidx=c(2001, 2006, 2011, 2016)- startyr+1L, agspan=5)
    fit$relincid <- estci2(abind::abind(lapply(fit$relincid, function(x) sweep(x, 2:3, x[3,,], "/")), rev.along=0))
    dimnames(fit$relincid)[1:3] <- list(agegr=paste0(3:9*5, "-", 3:9*5+4), sex=c("male", "female"), year=2001+0:3*5)
  }
  
  return(fit)
}

simfit.eppfit <- function(fit, rwproj=fit$fp$eppmod == "rspline", pregprev=TRUE){

  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    lastdata.idx <- as.integer(max(fit$likdat$ancsite.dat$df$yidx,
                                   fit$likdat$hhs.dat$yidx,
                                   fit$likdat$ancrtcens.dat$yidx,
                                   fit$likdat$hhsincid.dat$idx,
                                   fit$likdat$sibmx.dat$idx))

    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (lastdata.idx-1)/fit$fp$dt+1

    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- sim_rvec_rwproj(par$rvec, firstidx, lastidx, fit$fp$dt); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)
  
  fit$rvec <- sapply(mod.list, attr, "rvec")
  fit$prev <- sapply(mod.list, prev)
  fit$incid <- mapply(incid, mod = mod.list, fp = fp.list)
  fit$popsize <- sapply(mod.list, pop15to49)

  if(pregprev)
    fit$pregprev <- mapply(fnPregPrev, mod.list, fp.list)

  return(fit)
}


sim_mod_list <- function(fit, rwproj=fit$fp$eppmod == "rspline"){

  fit$param <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))

  if(rwproj){
    if(exists("eppmod", where=fit$fp) && fit$fp$eppmod == "rtrend")
      stop("Random-walk projection is only used with r-spline model")

    dt <- if(inherits(fit$fp, "eppfp")) fit$fp$dt else 1.0/fit$fp$ss$hiv_steps_per_year

    lastdata.idx <- as.integer(max(fit$likdat$ancsite.dat$df$yidx,
                                   fit$likdat$hhs.dat$yidx,
                                   fit$likdat$ancrtcens.dat$yidx,
                                   fit$likdat$hhsincid.dat$idx,
                                   fit$likdat$sibmx.dat$idx))
    
    fit$rvec.spline <- sapply(fit$param, "[[", "rvec")
    firstidx <- which(fit$fp$proj.steps == fit$fp$tsEpidemicStart)
    lastidx <- (lastdata.idx-1)/dt+1

    ## replace rvec with random-walk simulated rvec
    fit$param <- lapply(fit$param, function(par){par$rvec <- epp:::sim_rvec_rwproj(par$rvec, firstidx, lastidx, dt); par})
  }
  
  fp.list <- lapply(fit$param, function(par) update(fit$fp, list=par))
  mod.list <- lapply(fp.list, simmod)

  ## strip unneeded attributes to preserve memory

  keep <- c("class", "dim", "infections", "hivdeaths", "natdeaths", "hivpop", "artpop", "rvec", "popadjust")
  mod.list <- lapply(mod.list, function(mod){ attributes(mod)[!names(attributes(mod)) %in% keep] <- NULL; mod})

  return(mod.list)
}

## ' aggregate lists of model fits 
aggr_specfit <- function(fitlist, rwproj=sapply(fitlist, function(x) x$fp$eppmod) == "rspline"){
  allmod <- mapply(sim_mod_list, fitlist, rwproj, SIMPLIFY=FALSE)

  modaggr <- lapply(do.call(mapply, c(FUN=list, allmod, SIMPLIFY=FALSE)), Reduce, f="+")
  ##
  infectionsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "infections"), SIMPLIFY=FALSE)), Reduce, f="+")
  hivdeathsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "hivdeaths"), SIMPLIFY=FALSE)), Reduce, f="+")
  natdeathsaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "natdeaths"), SIMPLIFY=FALSE)), Reduce, f="+")
  hivpopaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "hivpop"), SIMPLIFY=FALSE)), Reduce, f="+")
  artpopaggr <- lapply(do.call(mapply, c(FUN=list, lapply(allmod, lapply, attr, "artpop"), SIMPLIFY=FALSE)), Reduce, f="+")
  ##
  modaggr <- mapply("attr<-", modaggr, "infections", infectionsaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "hivdeaths", hivdeathsaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "natdeaths", natdeathsaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "hivpop", hivpopaggr, SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "artpop", artpopaggr, SIMPLIFY=FALSE)
  ##
  modaggr <- mapply("attr<-", modaggr, "prev15to49", lapply(modaggr, calc_prev15to49, fitlist[[1]]$fp), SIMPLIFY=FALSE)
  modaggr <- mapply("attr<-", modaggr, "incid15to49", lapply(modaggr, calc_incid15to49, fitlist[[1]]$fp), SIMPLIFY=FALSE)
  ##
  modaggr <- lapply(modaggr, "class<-", c("specaggr", "spec"))
  return(modaggr)
}
