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

