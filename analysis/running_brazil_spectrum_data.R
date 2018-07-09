##############################################################################################################
## Running the model with Brazil data ########################################################################
##############################################################################################################

library(eppasm)
## devtools::load_all("~/Documents/Code/R/eppasm-csavr/") # @csavr
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

library(epp)

library(magrittr)
library(broom)
library(ggplot2)

##############################################################################################################
## Lets load up the data file and run the deterministic model ################################################
##############################################################################################################

brazil_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pjnz)
brazil_fp$artmx_timerr <- rep(1.0, brazil_fp$ss$PROJ_YEARS)
brazil_fp$t_diagn_start <- 16L   # assume diagnoses starts in 1985
brazil_fp$diagn_rate <- array(0.2, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_mod <- simmod(brazil_fp,VERSION = "R")
percent_undiag <- attr(brazil_mod,"undiagnosed_percent")
plot(percent_undiag,type="l",col="orange")

brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- F

#############################################################################################################
## Lets read in the csvar data for brazil from spectrum #####################################################
#############################################################################################################

brazil_csavrd <- read_csavr_data(brazil_pjnz)
brazil_csavrd[brazil_csavrd == 0] <- NA
brazil_csavrd$idx <- brazil_csavrd$year - brazil_fp$ss$proj_start + 1L

brazil <- list(fp = brazil_fp, csavrd = brazil_csavrd)
colnames(brazil$csavrd)[4] <- "total_cases"

brazil_opt1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, optfit=TRUE)
brazil_fit1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)


###########################################################################################################
## Alter the fp file to have four disease stages ##########################################################
###########################################################################################################

# 1 = 1 
# 2 = 2 
# 3 = 3 + 4
# 4 = 5 + 6 + 7

brazil_fp$ss$hDS <- 4L

brazil_fp$cd4_initdist[3,,] <- brazil_fp$cd4_initdist[3,,] + brazil_fp$cd4_initdist[4,,] 
brazil_fp$cd4_initdist[4,,] <- brazil_fp$cd4_initdist[5,,] + brazil_fp$cd4_initdist[6,,] + brazil_fp$cd4_initdist[7,,] 
brazil_fp$cd4_initdist <- brazil_fp$cd4_initdist[-c(5,6,7),,]

brazil_fp$cd4_prog[3,,] <- brazil_fp$cd4_prog[4,,]                                              ### Here the rate of leaving the 
brazil_fp$cd4_prog <- brazil_fp$cd4_prog[-c(4,5,6),,]                                           ### 3rd class is rate of leaving 
                                                                                                ### original's 4th
fp_alter <- function(fp_func,mean = F){
  if(nrow(fp_func) != 7 & nrow(fp_func) != 6 | length(dim(fp_func)) != 3){
    stop("Don't know how to deal with this input")
  }
  
  
  if(nrow(fp_func) == 7){
    if(mean == F){
    fp_func[3,,] <- fp_func[3,,] + fp_func[4,,]
    fp_func[4,,] <- fp_func[5,,] + fp_func[6,,] + fp_func[7,,]
    fp_func <- fp_func[-c(5,6,7),,]
    } else {
      fp_func[3,,] <- ((fp_func[3,,] * 2) + fp_func[4,,]) / 3
      fp_func[4,,] <- ((fp_func[5,,] * 2) + fp_func[6,,] + fp_func[7,,]) / 4
      fp_func <- fp_func[-c(5,6,7),,]
      
    }
  }
  if(nrow(fp_func) == 6){
    if(mean == F){
    fp_func[3,,] <- fp_func[3,,] + fp_func[4,,]
    fp_func <- fp_func[-c(4,5,6),,]
    }else{
      fp_func[3,,] <- (fp_func[3,,] + fp_func[4,,]) / 2
      fp_func <- fp_func[-c(4,5,6),,]
    }
  }
  return(fp_func)
}

brazil_fp$cd4_mort <- fp_alter(brazil_fp$cd4_mort, mean = T)
brazil_fp$frr_cd4 <- fp_alter(brazil_fp$frr_cd4)  
brazil_fp$paedsurv_cd4dist <- fp_alter(brazil_fp$paedsurv_cd4dist)

brazil_fp$frr_art[,3,,] <- brazil_fp$frr_art[,3,,] + brazil_fp$frr_art[,4,,]
brazil_fp$frr_art[,4,,] <- brazil_fp$frr_art[,5,,] + brazil_fp$frr_art[,6,,] + brazil_fp$frr_art[,7,,]
brazil_fp$frr_art <- brazil_fp$frr_art[,-c(5,6,7),,]

brazil_fp$art_mort[,3,,] <- (brazil_fp$art_mort[,3,,] * 2/3) + (brazil_fp$art_mort[,4,,] * 1/3)
brazil_fp$art_mort[,4,,] <- (brazil_fp$art_mort[,5,,] * 1/2) + (brazil_fp$art_mort[,6,,] * 1/4) + (brazil_fp$art_mort[,7,,] * 1/4) 
brazil_fp$art_mort <- brazil_fp$art_mort[,-c(5,6,7),,]

brazil_fp$paedsurv_artcd4dist[3,4,,] <- 1
brazil_fp$paedsurv_artcd4dist[1:2,,,] <- 0
brazil_fp$paedsurv_artcd4dist <- brazil_fp$paedsurv_artcd4dist[,-c(5,6,7),,]
brazil_fp$diagn_rate <- array(0.2, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))

brazil_fp$med_cd4init_cat 

brazil_mod <- simmod(brazil_fp,VERSION = "R")

percent_undiag <- attr(brazil_mod,"undiagnosed_percent")
lines(percent_undiag,col="midnightblue")

###################################################################################################
## Now we've modified the input data to have 4 cd4 classes let's use our Brazil dataset !!! #######
###################################################################################################

brazil_fp$time_at_which_get_cd4_counts <- 2001

brazil_deaths <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/Total_deaths_2017_2.csv")
brazil_detection <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/number_of_detected_hiv_cases_1980_2016_Brazil.csv")
brazil_cd4_at_detec <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/cd4_athiv_detection_brazil.csv")
brazil_art_number <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/siclom_number_of_new_ART_starters.csv")
brazil_cd4_art <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/CD4_at_ART_initiation.csv")

brazil_deaths_csav <- brazil_deaths[,c(1,8)]
brazil_deaths_csav$Corr_Total <- round(brazil_deaths_csav$Corr_Total)

brazil_csvard_test <- brazil_deaths_csav
names(brazil_csvard_test) <- c("year","aids_deaths")

detec_stage_1 <- round(brazil_cd4_at_detec$X.500[1:15] * brazil_detection$TOTAL[21:35])
detec_stage_2 <- round(brazil_cd4_at_detec$X350.500[1:15] * brazil_detection$TOTAL[21:35])
detec_stage_3 <- round(brazil_cd4_at_detec$X200.350[1:15] * brazil_detection$TOTAL[21:35])
detec_stage_4 <- round(brazil_cd4_at_detec$X.200[1:15] * brazil_detection$TOTAL[21:35])

brazil_csvard_test <- merge(brazil_csvard_test, brazil_detection[,c(1,4)], by = "year")
names(brazil_csvard_test) <- c("year","aids_deaths","total_cases")
brazil_csvard_test$stage_1_cases <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_1_cases[21:35] <- detec_stage_1

brazil_csvard_test$stage_2_cases <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_2_cases[21:35] <- detec_stage_2

brazil_csvard_test$stage_3_cases <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_3_cases[21:35] <- detec_stage_3

brazil_csvard_test$stage_4_cases <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_4_cases[21:35] <- detec_stage_4

names(brazil_art_number) <- c("year","total_art")
brazil_csvard_test$total_art <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$total_art[26:35] <- brazil_art_number$total_art[6:15]

art_stage_1 <- round(brazil_cd4_art$X.500[7:16] * brazil_art_number$total_art[6:15])
art_stage_2 <- round(brazil_cd4_art$X350.500[7:16] * brazil_art_number$total_art[6:15])
art_stage_3 <- round(brazil_cd4_art$X200.350[7:16] * brazil_art_number$total_art[6:15])
art_stage_4 <- round(brazil_cd4_art$X.200[7:16] * brazil_art_number$total_art[6:15])

brazil_csvard_test$stage_1_art <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_1_art[26:35] <- art_stage_1

brazil_csvard_test$stage_2_art <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_2_art[26:35] <- art_stage_2

brazil_csvard_test$stage_3_art <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_3_art[26:35] <- art_stage_3

brazil_csvard_test$stage_4_art <- rep(0, nrow(brazil_csvard_test))
brazil_csvard_test$stage_4_art[26:35] <- art_stage_4

brazil_csvard_test[brazil_csvard_test == 0] <- NA

brazil_csvard_test$idx <- brazil_csvard_test$year - brazil_fp$ss$proj_start + 1L


save(brazil_csvard_test,
     file="C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/my_created_csvard_for_brazil_WITH_ART")

###################################################################################################
## The moment of truth, let's try and run the brazil model with our csvar data ####################
###################################################################################################
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- F
brazil_fp$artcd4elig_idx <- rep(1L, length(brazil_fp$artcd4elig_idx))
brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

brazil_mod <- simmod(brazil$fp,VERSION = "R")
brazil_mod_out <- attributes(brazil_mod)
lines(brazil_mod_out$undiagnosed_percent,col="blueviolet")
mod_aids_deaths <- colSums(brazil_mod_out$hivdeaths,,2)


ll_aidsdeaths <- sum(dpois(brazil$csavrd$aids_deaths[4:35], mod_aids_deaths[15:46] , log=TRUE))
compo <- cbind.data.frame(brazil$csavrd$aids_deaths[4:35], mod_aids_deaths[15:46])
compo
mod_aid
brazil_opt1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, optfit=TRUE)
brazil_fit1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit double logistic model for incidence rate
brazil_opt2 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#############################################################################################################
## Now lets change it back to our 7 classes, the only slight of hand here is now the dfs for likelihood are #
## summed across the requisite classes ######################################################################
#############################################################################################################
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

brazil_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pjnz)
brazil_fp$artmx_timerr <- rep(1.0, brazil_fp$ss$PROJ_YEARS)
brazil_fp$t_diagn_start <- 10L   # assume diagnoses starts in 1985
brazil_fp$diagn_rate <- array(0.2, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_mod <- simmod(brazil_fp,VERSION = "R")
percent_undiag <- attr(brazil_mod,"undiagnosed_percent")
plot(percent_undiag,type="l",col="orange")

brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- F
brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- F
brazil_fp$time_at_which_get_cd4_counts <- 2001
brazil_fp$stages <- 4L
brazil_fp$artcd4elig_idx <- rep(1L, length(brazil_fp$artcd4elig_idx))
brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)
brazil_mod <- simmod(brazil$fp,VERSION = "R")
brazil_mod_out <- attributes(brazil_mod)
lines(brazil_mod_out$undiagnosed_percent,col="blueviolet")
mod_aids_deaths <- colSums(brazil_mod_out$hivdeaths,,2)


ll_aidsdeaths <- sum(dpois(brazil$csavrd$aids_deaths[4:35], mod_aids_deaths[15:46] , log=TRUE))
compo <- cbind.data.frame(brazil$csavrd$aids_deaths[4:35], mod_aids_deaths[15:46])
compo

brazil_opt1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit double logistic model for incidence rate
brazil_opt2 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

