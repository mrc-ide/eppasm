devtools::load_all("~/Dropbox/Documents/Code/R/eppspectrum/", export_all=FALSE)

upd.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Botswana_72.upd"
spec.path <- "~/Documents/Data/Spectrum files/2016 final/SSA/Botswana_ Final_15_04_ 2016 upd/Botswana_ Final_15_04_ 2016 upd"


demp <- read_specdp_demog_param(paste0(spec.path, ".DP"))
projp <- read_hivproj_param(paste0(spec.path, ".DP"))

fp <- create_spectrum_fixpar(projp, demp, proj_end = 2016, time_epi_start = projp$t0, hiv_steps_per_year= 5L)  # Set time_epi_start tomatch EPP


theta.rspline <- c(1.31820011, -0.09884313, -0.40054248, 0.06277183, 0.16923859, 0.41277390, -0.17640756, -14.13863910, 0.09765759, -3.73232668, -5.12046650)

param <- fnCreateParam(theta.rspline, fp)
fp <- update(fp, list=param)
fp$iota <- fp$iota*15

modR <- simmod(fp, VERSION="R")
modC <- simmod(fp)

round(prev(modR), 4)               # prevalence
round(prev(modC), 4)               # prevalence

round(incid(modR), 8)  # incidence
round(incid(modC), 8)  # incidence


library(microbenchmark)
mb <- microbenchmark

mb(simmod(fp), times=1e4)

## > mb(simmod(fp), times=1e4)
## Unit: microseconds
##        expr    min       lq     mean  median       uq      max neval
##  simmod(fp) 637.46 644.8285 779.8363 679.267 692.4155 52236.08 10000

## > mb(simmod(fp), times=1e4)
## Unit: microseconds
##        expr     min      lq     mean  median      uq      max neval
##  simmod(fp) 528.866 536.597 789.8745 594.865 746.289 50992.61 10000
