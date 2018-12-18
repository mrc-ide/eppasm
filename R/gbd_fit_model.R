## Prepare fit by EPP regions
#'
#' @param pjnz file path to Spectrum PJNZ file.
#' @param proj.end end year for projection.
#' @param popupdate logical should target population be updated to match
#'   age-specific population size from DP file and %Urban from EPP XML.
#'   
prepare_spec_fit_gbd <- function(loc, collapse, i, proj.end=2016.5, gbd.pop = TRUE, popadjust = NULL, popupdate=TRUE, use_ep5=FALSE){
  
  pjnz <- find_pjnz(loc)[[1]]
  ## epp
  if(!collapse | grepl ('ZAF', loc)){
  eppd <- epp::read_epp_data(pjnz)
  epp.subp <- epp::read_epp_subpops(pjnz)
  epp.input <- epp::read_epp_input(pjnz)
    if(grepl('ZAF', loc)){
      zaf.dict <- list("MP" = "ZAF_487", "GP" = "ZAF_484", "KZN" = "ZAF_485", 
                       "WC" = "ZAF_490", "EC" = "ZAF_482", "LP" = "ZAF_486", 
                       "FS" = "ZAF_483", "NW" = "ZAF_488", "NC" = "ZAF_489")
      eppd.new <- list()
      eppd.new[[loc]] <- eppd[[names(which(zaf.dict == loc))]]  
      eppd <- eppd.new
      epp.subp$total <- epp.subp$subpops[[names(which(zaf.dict == loc))]] 
      eppd[[1]]$ancrtcens <- data.frame(year=integer(), prev=integer(), n=integer())
        # eppd.tot[[subpop.tot]]$ancrtcens <- NULL
      epp.subp$subpops <- NULL
      epp.subp$subpops[[loc]] <- epp.subp$total
      ## TODO: What to do with epp.input? (especially ART values)
    }
  } else{
    epp.totals <- collapse_epp(loc)
    eppd <- epp.totals$eppd
    epp.subp <- epp.totals$epp.subp.tot
    epp.input <- epp.totals$epp.input.tot
  }
  if(gbd.pop){
    epp.pops <- sub.pop.params.epp(epp.subp, epp.input, loc)
    epp.subp <- epp.pops[['epp.subp']]
    epp.input <- epp.pops[['epp.input']]
  }

  epp.subp.input <- epp::fnCreateEPPSubpops(epp.input, epp.subp, eppd)

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
  if(gbd.pop){
    demp <- sub.pop.params.demp(demp, loc, i)
  }
  projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)
  # epp_t0 <- read_epp_t0(pjnz)
  ## TODO: figure out what to do with t0
  epp_t0 <- epp.input$epidemic.start
  ## If popadjust = NULL, look for subp if more than 1 EPP region
  if(is.null(popadjust))
    popadjust <- length(eppd) > 1

  ## If Urban/Rural fit, read percentage urban from EPP XML file
  if(length(eppd) == 2 && all(sort(substr(names(eppd), 1, 1)) == c("R", "U")))
    perc_urban <- read_epp_perc_urban(pjnz) / 100
  else
    perc_urban <- NULL
    
  ##TODO: don't need this if running for gbd pop?
  specfp.subp <- create_subpop_specfp(projp, demp, eppd, proj_end=proj.end, epp_t0=epp_t0,
                                      popadjust = popadjust, popupdate = popupdate, perc_urban = perc_urban)
  
  if(is.na(specfp.subp[[loc]]$ss$time_epi_start)) {specfp.subp[[loc]]$ss$time_epi_start <- epp_t0}
  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, epp::fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "specfp", specfp.subp)
  val <- set.list.attr(val, "country", read_country(pjnz))
  val <- set.list.attr(val, "region", names(eppd))

  attr(val, "country") <- read_country(pjnz)
  attr(val, "region") <- read_region(pjnz)

  return(val)
}

