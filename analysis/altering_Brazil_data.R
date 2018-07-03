#####################################################################################################################################
## Adapting chile to Brazil population size, very dodgy!! ###########################################################################
#####################################################################################################################################

library(knitr)
opts_chunk$set(tidy=TRUE, warning=FALSE, cache=TRUE, message=FALSE)
options(knitr.kable.NA = '')
library(eppasm)
## devtools::load_all("~/Documents/Code/R/eppasm-csavr/") # @csavr
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

library(epp)

library(magrittr)
library(broom)
library(ggplot2)

######################################################################################
## Now lets load up the chile .pngz file #############################################
######################################################################################

cl_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Chile_2017_final.pjnz"

cl_fp <- prepare_directincid(cl_pjnz)
cl_fp$artmx_timerr <- rep(1.0, cl_fp$ss$PROJ_YEARS)
cl_fp$t_diagn_start <- 16L   # assume diagnoses starts in 1985
cl_fp$diagn_rate <- array(0.2, c(dim(nl_fp$cd4_mort), cl_fp$ss$PROJ_YEARS))
cl_fp$relinfectART <- 0.3
cl_fp$tsEpidemicStart <- 1970.5

cl_mod <- simmod(cl_fp)

pop_difference <- 207.7 / 17.91

brazil_fp <- cl_fp
brazil_fp$basepop <- brazil_fp$basepop * pop_difference 

brazil_fp$birthslag <- brazil_fp$birthslag * pop_difference

brazil_fp$entrantpop <- brazil_fp$entrantpop * pop_difference
 
brazil_fp$targetpop <- brazil_fp$targetpop * pop_difference

brazil_fp$births <- brazil_fp$births * pop_difference

brazil_fp$art15plus_num <- brazil_fp$art15plus_num * pop_difference

###### So that's all the absolute numbers inflated with Brazil's increased population size 

brazil_mod <- simmod.specfp(brazil_fp,VERSION = "R")
percent_undiag_br <- attr(brazil_mod,"undiagnosed")

plot(percent_undiag_br,type="l",col="forestgreen")
abline(v=15,col="red")
abline(v=25,col="red")

#######################################################################################
###### Now lets adjust the casavrd data ###############################################
#######################################################################################

cl_csavrd <- read_csavr_data(cl_pjnz)
cl_csavrd[cl_csavrd == 0] <- NA
cl_csavrd$idx <- cl_csavrd$year - cl_fp$ss$proj_start + 1L

cl_csavrd
brazil_csvard <- cl_csavrd

brazil_cases <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/AIDS_cases_from_SINAN.csv")
brazil_cases <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/number_of_detected_hiv_cases_1980_2016_Brazil.csv")

for(i in 1:nrow(brazil_csvard)){
if( brazil_csvard[i,1] %in% brazil_cases[,1] == T){
  year <- brazil_csvard[i,1]
  brazil_csvard[i,4] <- brazil_cases[brazil_cases$year==year,3]
}
}

brazil_csvard
brazil_csvard[,-c(4)] <- NA

################################################################################################
####### Now we only have number of new cases, so lets just use that and see what happens ... ###
################################################################################################

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard)
brazil_opt1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, optfit=TRUE)
brazil_opt2 <- fitmod_csavr(cl, incid_func = "idbllogistic", B0=1e3, optfit=TRUE)

#################################################################################################
## Adjusting the likelihood for diagnoses to be by cd4 stage ####################################
#################################################################################################

brazil_mod_out <- attributes(brazil_mod)

diagnoses_rates_by_cd4_stage <- function(mod){

  diagnoses_full <- attributes(mod)
  
  diagnoses <- diagnoses_full$diagnoses

  disease_stage <- as.character(1:nrow(diagnoses))

  diagnoses_df <- NULL
for(i in 1:nrow(diagnoses)){
  mean_per_stage <- NULL
  for(ii in 1:52){
    mean_per_stage[ii] <- mean(diagnoses[i,,,ii])
  }
  mean_per_stage <- cbind.data.frame(mean_per_stage,c(1970:2021),rep(disease_stage[i],length(mean_per_stage)))
  diagnoses_df <- rbind.data.frame(diagnoses_df,mean_per_stage)
}
  names(diagnoses_df) <- c("diagnoses","year","stage")
  return(diagnoses_df)
}
  
diagnoses_by_stage <- diagnoses_rates_by_cd4_stage(brazil_mod)  

cl_csavrd

cl_mod_out <- attributes(cl_mod)
cl_mod_out

brazil_mod_out$diagnoses

mod_tot <- colSums(cl_mod_out$diagnoses,,3)

mod_tot_tims <- cbind.data.frame(mod_tot,1970:2021)
mod_tot_tims$cl_dat <- rep(0, nrow(mod_tot_tims))

for(i in 1:nrow(mod_tot_tims)){
  
  if(mod_tot_tims[i,2] %in% cl_csavrd$year){
    year <- mod_tot_tims[i,2]
    mod_tot_tims[i,3] <- cl_csavrd[cl_csavrd$year == year, 4] 
    
    print(year)
    Sys.sleep(0.5)
  }
  
}
mod_tot_tims
brazil_cases[brazil_cases$year == 1985,4]


sum(dpois(mod_tot_tims$cl_dat[1:47], round(mod_tot_tims$mod_tot[1:47]), log=TRUE))

##########################################################################################
## Working out how to fit with ART initiation data as well ###############################
##########################################################################################

dim(brazil_mod_out$artinits)

art_starts_mod <- colSums(brazil_mod_out$artinits,,3)

art_starts_to_fit_with <- art_starts_mod[(2006-1970):(2015-1970)]

brazil_art_starts <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/siclom_number_of_new_ART_starters.csv")
brazil_art_starts_to_fit_with <- brazil_art_starts[6:15,2]

sum(dpois(brazil_art_starts_to_fit_with,art_starts_to_fit_with,log = T))

mod_tot_tims$braz_dat[1:47]
round(mod_tot_tims$mod_tot[1:47])
