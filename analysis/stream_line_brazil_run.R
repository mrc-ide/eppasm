############################################################################################################
## Running brazil data 2 ###################################################################################
############################################################################################################

library(eppasm)
## devtools::load_all("~/Documents/Code/R/eppasm-csavr/") # @csavr
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

library(epp)

library(magrittr)
library(broom)
library(ggplot2)

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/my_created_csvard_for_brazil_WITH_ART",verbose = T)

############################################################################################################
## NOw lets run the thing !!! ##############################################################################
############################################################################################################

brazil_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pjnz)
brazil_fp$artmx_timerr <- rep(1.0, brazil_fp$ss$PROJ_YEARS)
brazil_fp$t_diagn_start <- 10L   # assume diagnoses starts in 1985


brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- F
brazil_fp$time_at_which_get_cd4_counts <- 2001

brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$likelihood_cd4 <- T
brazil_fp$artinit_use <- F
brazil_fp$artcd4elig_idx <- rep(1L, length(brazil_fp$artcd4elig_idx))
brazil_fp$stages <- 4L

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

#########################################################################################################
## Now lets get the right ART data in there #############################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/corrected_art_starters_BRAZIL",verbose = T)

brazil$fp$art15plus_num <- art_format

#########################################################################################################
## Now lets run the model with the raw diagnoses and AIDS deaths and see the output #####################
#########################################################################################################
brazil$fp$aidsdeaths <- T
brazil$fp$diagnoses_uses <- T
brazil$fp$tARTstart <- 28L
brazil$fp$linear_diagnosis <- TRUE

brazil_opt1_cd4_both <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1_cd4_both <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_cd4_both <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2_cd4_both <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_cd4_both <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3_cd4_both <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

save(brazil_fit1_cd4_both,file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/brazil_fit1_cd4_both")
save(brazil_fit2_cd4_both, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/brazil_fit2_cd4_both")
#### lets tidy the results and plot them #####
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/brazil_fit1_cd4_both",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/brazil_fit2_cd4_both", verbose = T)

brazil_out1_cd4_both <- tidy(brazil_fit1_cd4_both) %>% data.frame(model = "logistic", .)
brazil_out2_cd4_both <- tidy(brazil_fit2_cd4_both) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4_both <- tidy(brazil_fit3_cd4_both) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4_both, brazil_out2_cd4_both, brazil_out3_cd4_both)


no_cd4_fit <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=1)

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/CD4_both_DIAG_DEATHS_tara_ART")
save(brazil_out_cd4, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Data/CD4_both_DIAG_DEATHS_tara_ART_DATA")

########################################################################################################
## NOw lets just fit it with Diagnoses #################################################################
########################################################################################################
brazil$fp$aidsdeaths <- T
brazil$fp$diagnoses_uses <- F

brazil_opt1_no_diag <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1__no_diag <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_no_diag <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2_no_diag <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_no_diag <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3_no_diag <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#### lets tidy the results and plot them #####

brazil_out1_no_diag <- tidy(brazil_fit1__no_diag) %>% data.frame(model = "logistic", .)
brazil_out2_no_diag <- tidy(brazil_fit2_no_diag) %>% data.frame(model = "double logistic", .)
brazil_out3_no_diag <- tidy(brazil_fit3_no_diag) %>% data.frame(model = "rlogistic", .)
brazil_out_no_diag <- rbind(brazil_out1_no_diag, brazil_out2_no_diag, brazil_out3_no_diag)

brazil_out_no_diag[is.infinite(brazil_out_no_diag$mean)] <- 0
brazil_out_no_diag[brazil_out_no_diag$model =="logistic",]

brazil_out_no_diag_no_logi <- rbind(brazil_out2_no_diag, brazil_out3_no_diag)

no_cd4_fit <- ggplot(subset(brazil_out_no_diag, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  scale_y_continuous(breaks = NULL)

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_diag_tara_art_no_cd4_plot")
save(brazil_out_no_diag,file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Data/no_diag_tara_art_no_cd4_data")

#######################################################################################################
## Now we will run the data with just the diagnoses ###################################################
#######################################################################################################
brazil$fp$aidsdeaths <- F
brazil$fp$diagnoses_uses <- T

brazil_opt1_no_aids <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1__no_aids <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_no_aids <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2_no_aids <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_no_aids <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3_no_aids <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#### lets tidy the results and plot them #####

brazil_out1_no_aids <- tidy(brazil_fit1__no_aids) %>% data.frame(model = "logistic", .)
brazil_out2_no_aids <- tidy(brazil_fit2_no_aids) %>% data.frame(model = "double logistic", .)
brazil_out3_no_aids <- tidy(brazil_fit3_no_aids) %>% data.frame(model = "rlogistic", .)
brazil_out_no_aids <- rbind(brazil_out1_no_aids, brazil_out2_no_aids, brazil_out3_no_aids)

no_cd4_fit <- ggplot(subset(brazil_out_no_aids, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) #+
  # scale_y_continuous(breaks = NULL)

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_tara_art_no_cd4_plot")
save(brazil_out_no_diag,file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Data/no_AIDSDEATHS_tara_art_no_cd4_data")

#####################################################################################################################
## Now lets run it with CD4 counts for diagnoses only as well #######################################################
#####################################################################################################################

brazil$fp$likelihood_cd4 <- T

brazil_opt1_no_aids_cd4 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1__no_aids_cd4 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_no_aids_cd4 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2_no_aids_cd4 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_no_aids_cd4 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3_no_aids_cd4 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#### lets tidy the results and plot them #####

brazil_out1_no_aids_cd4 <- tidy(brazil_fit1__no_aids_cd4) %>% data.frame(model = "logistic", .)
brazil_out2_no_aids_cd4 <- tidy(brazil_fit2_no_aids_cd4) %>% data.frame(model = "double logistic", .)
brazil_out3_no_aids_cd4 <- tidy(brazil_fit3_no_aids_cd4) %>% data.frame(model = "rlogistic", .)
brazil_out_no_aids_cd4 <- rbind(brazil_out1_no_aids_cd4, brazil_out2_no_aids_cd4, brazil_out3_no_aids_cd4)

no_cd4_fit <- ggplot(subset(brazil_out_no_aids_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=1) #+
# scale_y_continuous(breaks = NULL)

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_CD4_likel_tara_art_no_cd4_plot")
save(brazil_out_no_diag,file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Data/no_AIDSDEATHS_CD4_likek_tara_art_no_cd4_data")

########################################################################################
## load up the plots to save ###########################################################
########################################################################################

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_CD4_likel_tara_art_no_cd4_plot",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_tara_art_no_cd4_plot", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_diag_tara_art_no_cd4_plot", verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/",verbose = T)
no_cd4_fit <- ggplot(subset(no_cd4_fit$data, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=1) +
  scale_y_continuous(breaks=NULL)
