prepare_spec_object <- function(loc, popadjust = TRUE, popupdate=TRUE, use_ep5=FALSE){
  
  pjnz <- find_pjnz(loc)[[1]]
  ##TODO: Make this work for ZAF
  if(grepl ('ZAF', loc)){
    eppd <- epp::read_epp_data(pjnz)
    zaf.dict <- list("MP" = "ZAF_487", "GP" = "ZAF_484", "KZN" = "ZAF_485", 
                     "WC" = "ZAF_490", "EC" = "ZAF_482", "LP" = "ZAF_486", 
                     "FS" = "ZAF_483", "NW" = "ZAF_488", "NC" = "ZAF_489")
    eppd.new <- list()
    eppd.new[[loc]] <- eppd[[names(which(zaf.dict == loc))]]  
    eppd <- eppd.new
    # eppd[[1]]$ancrtcens <- data.frame(year=integer(), prev=integer(), n=integer())
    # eppd.tot[[subpop.tot]]$ancrtcens <- NULL
  } else{
    epp.totals <- collapse_epp(loc)
    eppd <- epp.totals$eppd
  }
  
  country <- attr(eppd, "country")
  cc <- attr(eppd, "country_code")
  
  ## melt site-level data
  eppd <- Map("[[<-", eppd, "ancsitedat", lapply(eppd, melt_ancsite_data))
  
  ## tidy HHS data
  eppd <- Map("[[<-", eppd, "hhs", lapply(eppd, tidy_hhs_data))
  
  attr(eppd, "country") <- country
  attr(eppd, "country_code") <- cc
  eppd <- eppd[[loc]]
  
  ## spectrum
  ## TODO: Do we want to implement a collapse function for demog param? (probably not, but it's worth noting
  ## that for locations we collapse (like Benin, Cote d'Ivoire, etc), the demp object will only hold the population for 1 subnational)
  demp <- read_specdp_demog_param(pjnz, use_ep5=use_ep5)

  projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)
  epp_t0 <- epp.totals$epp.input.tot$epidemic.start

  ## If popadjust = NULL, look for subp if more than 1 EPP region
  if(is.null(popadjust)){
    popadjust <- length(eppd) > 1
  }

  
  ## TODO: Run this after subbing into projp and demp
  specfp <- create_spectrum_fixpar(projp, demp, popadjust=popadjust, time_epi_start=epp_t0)
  
  specfp$ss$time_epi_start <- epp_t0
  ## output
  val <- list()
  attr(val, 'eppd') <- eppd
  attr(val, 'specfp') <- specfp
  attr(val, 'country') <- read_country(pjnz)
  attr(val, 'region') <- loc
  saveRDS(val, paste0('/share/hiv/data/PJNZ_EPPASM_prepped/', loc, '.rds'))
  
  return(val)
}