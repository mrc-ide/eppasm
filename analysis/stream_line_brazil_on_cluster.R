#############################################################################################################
## Streamline cluster running ###############################################################################
#############################################################################################################

setwd("X:")
options(didehpc.username = "jd2117",didehpc.home = "X:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
didehpc::didehpc_config()

context::context_log_start()

root <- "contexts_2"

ctx<- context::context_save(root, packages = c("epp","eppasm","anclik"), 
                            package_sources = provisionr::package_sources(local = c("C:/Users/josh/Dropbox/hiv_project/eppasm",
                                                                          "C:/Users/josh/Documents/R/win-library/3.4/epp",
                                                                          "C:/Users/josh/Documents/R/win-library/3.4/anclik")))
config <- didehpc::didehpc_config(cores = 3, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)

##########################################################################################################################
### NOw lets load up the data to put into the cluster run ################################################################
##########################################################################################################################

obj$enqueue(devtools::load_all("X:/eppasm"))
obj$enqueue(devtools::build("X:/eppasm"))

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
brazil_fp$likelihood_cd4 <- F
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
brazil_fp$tARTstart <- 28


brazil_opt1_cd4 <-  obj$enqueue(eppasm::fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),name = "opt_fit_1_")
brazil_opt1_cd4$status()

brazil_fit1_cd4 <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_cd4 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)
brazil_fit2_cd4 <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_cd4 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)
brazil_fit3_cd4 <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#### lets tidy the results and plot them #####

brazil_out1_cd4 <- tidy(brazil_fit1_cd4) %>% data.frame(model = "logistic", .)
brazil_out2_cd4 <- tidy(brazil_fit2_cd4) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4 <- tidy(brazil_fit3_cd4) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4)



no_cd4_fit <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5)

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/PLOT_no_cd4_tara_ART_both_DIAG_and_AIDSdeaths")
save(brazil_out_cd4, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/DATA_no_cd4_tara_ART_both_DIAG_and_AIDSdeaths")



