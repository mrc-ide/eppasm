read_spec_object <- function(loc, i, start.year = 1970, stop.year = 2019, trans.params.sub = TRUE, 
                             pop.sub = TRUE, anc.sub = FALSE, prev.sub = TRUE, popadjust = TRUE, age.prev = FALSE){
  dt <- readRDS(paste0('/share/hiv/data/PJNZ_EPPASM_prepped/', loc, '.rds'))
  
  ## Substitute IHME data
  ## Population parameters
  if(pop.sub){
    ## Unfortunately currently necessary to read PJNZ here to sub directly into demp and then call create_spectrum_fixpar
    pjnz <- find_pjnz(loc)[[1]]
    demp <- read_specdp_demog_param(pjnz, use_ep5=use_ep5)
    demp <- sub.pop.params.demp(demp, loc, i)
    projp <- read_hivproj_param(pjnz, use_ep5=use_ep5)
    specfp <- create_spectrum_fixpar(projp, demp, proj_start = start.year, proj_end = stop.year, popadjust=popadjust, time_epi_start=attr(dt, 'specfp')$ss$time_epi_start)
    attr(dt, 'specfp') <- specfp
  }
  
  ## Prevalence surveys
  if(prev.sub) {
      print("Substituting prevalence surveys")
    if(age.prev){
      dt <- sub.prev.granular(dt, loc)
    } else{
      dt <- sub.prev(loc, dt)	
    }
  }
  
  ## ANC data
  # high.risk.list <- loc.table[epp == 1 & collapse_subpop == 0 & !grepl("ZAF", ihme_loc_id) & !grepl("KEN", ihme_loc_id), ihme_loc_id]
  # ken.anc.path <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/data/prepped/kenya_anc_map.csv")
  # ken.anc <- fread(ken.anc.path)
  # no.anc.ken <- setdiff(loc.table[epp == 1 & grepl("KEN", ihme_loc_id), ihme_loc_id], ken.anc$ihme_loc_id)
  # if(loc %in% c(high.risk.list, "PNG", no.anc.ken)) {
  # 	anc.backcast <- F
  # }
  # if(anc.backcast) {
  # 	dt <- sub.anc(loc, dt)
  # }
  
  ## Transition parameters
  if(trans.params.sub) {
    print('Substituting transition parameters')
    dt <- sub.off.art(dt, loc, i)
    dt <- sub.on.art(dt, loc, i)
    dt <- sub.cd4.prog(dt, loc, i)
  }
  
}
  