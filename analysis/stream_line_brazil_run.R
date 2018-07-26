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
brazil_fp$likelihood_cd4 <- FALSE
brazil_fp$artinit_use <- FALSE
brazil_fp$time_at_which_get_cd4_counts <- 2001

brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$likelihood_cd4 <- F
brazil_fp$artinit_use <- FALSE
brazil_fp$artcd4elig_idx <- rep(1L, length(brazil_fp$artcd4elig_idx))
brazil_fp$stages <- 4L

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

#########################################################################################################
## Now lets get the right ART data in there #############################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/corrected_art_starters_BRAZIL",verbose = T)

brazil$fp$art15plus_num <- art_format

#########################################################################################################
## Lets use tara's art mortality rates ##################################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/corrected_ART_mort_brazil", verbose = T)

brazil$fp$art_mort <- test_art_mort

#########################################################################################################
## Now lets run the model with the raw diagnoses and AIDS deaths and see the output #####################
#########################################################################################################
brazil$fp$aidsdeaths <- FALSE
brazil$fp$diagnoses_uses <- TRUE
brazil$fp$tARTstart <- 28L
brazil$fp$likelihood_cd4 <- TRUE
brazil$fp$linear_diagnosis <- "12 param"




brazil_opt1_cd4_both <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)
brazil_fit1_cd4_both <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_opt2_cd4_both <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e3, optfit=TRUE)
brazil_fit2_cd4_both <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
brazil_opt3_cd4_both <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e3, optfit=TRUE)
brazil_fit3_cd4_both <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

brazil_out1_cd4_both <- tidy(brazil_fit1_cd4_both) %>% data.frame(model = "logistic", .)
brazil_out2_cd4_both <- tidy(brazil_fit2_cd4_both) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4_both <- tidy(brazil_fit3_cd4_both) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4_both, brazil_out2_cd4_both, brazil_out3_cd4_both)


fit_output <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=1)

save(fit_output, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/CD4_both_DIAG_DEATHS_tara_ART")
save(brazil_out_cd4, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Data/CD4_both_DIAG_DEATHS_tara_ART_DATA")

########################################################################################################
## NOw lets just fit it with Diagnoses #################################################################
########################################################################################################
brazil$fp$aidsdeaths <- TRUE
brazil$fp$diagnoses_uses <- FALSE

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


no_diag_fit <- ggplot(subset(brazil_out_no_diag, year %in% 1975:2017), 
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
brazil$fp$aidsdeaths <- FALSE
brazil$fp$diagnoses_uses <- TRUE

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

brazil$fp$likelihood_cd4 <- TRUE

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


##############################################################################
## Optimizer fiddle ##########################################################
##############################################################################

brazil_opt_2_out <- attributes(brazil_opt2_cd4_both$mod)
hiv_deaths <- colSums(brazil_opt_2_out$hivdeaths,,2)
plot(hiv_deaths,type="l")
points(brazil_opt2_cd4_both$dat$aids_deaths,col="red")

diagnoses <- colSums(brazil_opt_2_out$diagnoses,,3)
plot(diagnoses, type = "l", col="blue")
indexed_diag <- c(rep(0,10),brazil_opt2_cd4_both$dat$total_cases[1],0,brazil_opt2_cd4_both$dat$total_cases[2:35])
points(indexed_diag,col="red")


##############################################################################
## loading up the plots ######################################################
##############################################################################

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/CD4_both_DIAG_DEATHS_tara_ART",verbose = T)
CD4_both_diag_and_Deaths_fit <- no_cd4_fit
CD4_both_diag_and_Deaths_fit

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_CD4_likel_tara_art_no_cd4_plot", verbose = T)
cd4_diag_no_deaths_fit <- no_cd4_fit
cd4_diag_no_deaths_fit

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_AIDSDEATHS_tara_art_no_cd4_plot", verbose = T)
no_cd4_diag_no_deaths_fit <- no_cd4_fit
no_cd4_diag_no_deaths_fit

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/no_diag_tara_art_no_cd4_plot", verbose = T)
no_diag_DEATHS_fit <- no_cd4_fit
no_diag_DEATHS_fit

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/Plots/PLOT_no_cd4_tara_ART_both_DIAG_and_AIDSdeaths", verbose = T)
no_cd4_diag_and_deaths_fit <- no_cd4_fit
no_cd4_diag_and_deaths_fit


############################################################################################################
## Runnning the optims for the aidsdeaths, lets use the 6 knot linear diagnosis rate #######################
############################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"

brazil_opt1_knot_linear <- fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE)

brazil_opt2_knot_linear <- fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE)

brazil_opt3_knot_linear <- fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE)

###########################################################################################################
## Runnning the debugger on the eppasm R framework to see if we can hard code the numbers initating ART ###
## by CD4 count from the data that we have ################################################################
###########################################################################################################
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")

load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_knotlinear_mod_2",verbose = T)
load("C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_knotlinear_mod_3", verbose = T)

opt_2_res$fp
opt_3_res$fp


tara_art_cd4 <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/CD4_at_ART_initiation.csv")

art_cd4_counts <- t(tara_art_cd4)
art_cd4_counts <- art_cd4_counts[-1,-17]
rownames(art_cd4_counts) <- c()
class(art_cd4_counts)

save(art_cd4_counts, 
     file = "C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/cd4_at_ART_MATRIX")

opt_2_res$fp$tara_art_cd4 <- art_cd4_counts
opt_2_res$fp$use_cd4_art <- TRUE
opt_2_res$fp$time_cd4_stop <- 2015L
opt_2_res$fp$artcd4elig_idx <- rep(1L, 52)

opt_3_res$fp$tara_art_cd4 <- art_cd4_counts
opt_3_res$fp$use_cd4_art <- TRUE
opt_3_res$fp$time_cd4_stop <- 2015L
opt_3_res$fp$artcd4elig_idx <- rep(1L, 52)

test_run_opt_1 <- simmod.specfp(opt_2_res$fp, VERSION = "R")

brazil$fp$incidinput <- opt_3_res$fp$incidinput
brazil$fp$diagn_rate <- opt_3_res$fp$diagn_rate
brazil$fp$tara_art_cd4 <- art_cd4_counts
brazil$fp$use_cd4_art <- TRUE
brazil$fp$time_cd4_stop <- 2015L
brazil$fp$artcd4elig_idx <- rep(1L, 52)

test_run_opt_2 <- simmod.specfp(brazil$fp, VERSION = "R")

opt_2_res$fp$tara_art_cd4 <- art_cd4_counts
opt_2_res$fp$use_cd4_art <- FALSE
opt_2_res$fp$time_cd4_stop <- 2015L
opt_2_res$fp$artcd4elig_idx <- rep(1L, 52)

test_run_opt_3 <- simmod.specfp(opt_2_res$fp, VERSION = "R")

opt_2_res$fp$use_cd4_art <- FALSE
opt_2_res$fp$cd4_mort <- opt_3_res$fp$cd4_mort 
opt_2_res$fp$cd4_prog <- opt_3_res$fp$cd4_prog 
opt_2_res$fp$artcd4elig_idx <- rep(7, 52)



test_run_opt_4 <- simmod.specfp(opt_2_res$fp, VERSION = "R")

plot_undiagnosed <- function(optim_output,diag_start = 1980, art_start = 1996, mod_output = FALSE){
  tot_undiag_data <- NULL
  input_label <- paste("input",(1:length(optim_output)),sep = " ")
  tot_deaths <- NULL
  tot_diagnoses <- NULL
  tot_diag_rate <- NULL
  art_init_tot <- NULL
  tot_divided <- NULL
  
  
  for(i in 1:length(optim_output)){
    if(length(optim_output) == 1){
      if(mod_output == FALSE){
      list_version <- attributes(optim_output[[1]]$mod)
      }else{
        list_version <- attributes(optim_output)
      }
    }else{
      if(mod_output == FALSE){
      list_version <- attributes(optim_output[[i]]$mod)
      }else
        list_version <- attributes(optim_output[[i]])
    }
    undiag_percent <- (apply(list_version$hivpop,4,sum) - apply(list_version$diagnpop,4,sum)) / 
      (apply(list_version$hivpop,4,sum) + apply(list_version$artpop,5,sum)) * 100
    
    percent_dat <- cbind.data.frame(undiag_percent,rep(input_label[i],length(undiag_percent)),c(1970:2021))
    ## deaths
    
    deaths <- colSums(list_version$hivdeaths,,2)
    deaths_df <- cbind.data.frame(deaths, rep(input_label[i], length(deaths)), c(1970:2021))
    
    tot_deaths <- rbind.data.frame(tot_deaths, deaths_df)
    ## diagnoses
    
    diagnoses <- colSums(list_version$diagnoses,, 3)
    diagnoses_df <- cbind.data.frame(diagnoses, rep(input_label[i], length(diagnoses)), c(1970:2021))
    
    tot_diagnoses <- rbind.data.frame(tot_diagnoses, diagnoses_df)
    
    tot_undiag_data <- rbind.data.frame(tot_undiag_data,percent_dat)
    
    ## Diag rate 
    if(mod_output == FALSE){
    diag_rate <- optim_output[[i]]$par[(length(optim_output[[i]]$par) - 4) : length(optim_output[[i]]$par)]
    knots <- c(1986, 1996, 2000, 2009, 2015)
    
    diagn_trend <- approx(knots, diag_rate, 1970:2021, rule = 2)$y
    diagn_df <- cbind.data.frame(diagn_trend, rep(input_label[i], 52), c(1970:2021))
    
    tot_diag_rate <- rbind.data.frame(tot_diag_rate, diagn_df)
    }
    ## ARt init cats
    get_dist_per_year <- array(0, dim = c(7,52))
    get_dist_4_cat <- array(0, dim = c(4, 52))
    
    for(j in 1:52){
      for(ii in 1:7){
        
        get_dist_per_year[ii, j] <- sum(list_version$artinits[ii,,,j])
        
      }
      get_dist_4_cat[1,j] <- get_dist_per_year[1,j]
      get_dist_4_cat[2, j] <- get_dist_per_year[2, j]
      get_dist_4_cat[3, j] <- get_dist_per_year[3, j] + get_dist_per_year[4, j]
      get_dist_4_cat[4, j] <- get_dist_per_year[5, j] + get_dist_per_year[6, j] + get_dist_per_year[7, j]
      
      
      
    }
    
    get_dist_4_cat <- t(get_dist_4_cat)
    get_dist_4_cat <- data.frame(get_dist_4_cat)
    get_dist_4_cat$time <- c(1970:2021)
    
    funky_dist_4 <- reshape2::melt(get_dist_4_cat, id = "time")
    
    funky_dist_4$input <- rep(input_label[i], nrow(funky_dist_4))
    
    
    art_init_tot <- rbind.data.frame(art_init_tot, funky_dist_4)
    
    ### Deaths by wether on ART on not
    
    deaths_by_hiv <- colSums(list_version$hivpopdeaths,,3)
    deaths_on_art <- colSums(list_version$artpopdeaths,,4)
    
    deaths_divided_df <- cbind.data.frame(deaths_by_hiv, deaths_on_art, c(1970:2021))
    names(deaths_divided_df) <- c("hiv","art","time")
    
    deaths_funky <- reshape2::melt(deaths_divided_df, id = "time")
    deaths_funky$input <- rep(input_label[i], nrow(deaths_funky))
    
    tot_divided <- rbind.data.frame(tot_divided, deaths_funky)
    
    
    
    
    
  }
  
  names(tot_undiag_data) <- c("undiagnosed","input","year")
  names(tot_deaths) <- c("deaths", "input", "year")
  names(tot_diagnoses) <-c("diagnoses", "input", "year")
  if(mod_output == FALSE){
  names(tot_diag_rate) <-c("rate","input","year")
  }
  names(art_init_tot) <- c("time","cat","val","input")
  names(tot_divided) <- c("time","class","deaths","input")
  
  art_init_tot$cat_input <- paste(art_init_tot$cat, art_init_tot$input, sep = "")
  tot_divided$class_input <- paste(tot_divided$class, tot_divided$input, sep ="")
  
  
  
  if(mod_output == FALSE){
  tot_deaths <- merge(tot_deaths, optim_output[[1]]$likdat$aidsdeaths, all = TRUE)
  tot_diagnoses <- merge(tot_diagnoses, optim_output[[1]]$likdat$diagnoses, all = TRUE)
  }
  
  undiag_plot <- ggplot(data = tot_undiag_data,aes(x=year,y=undiagnosed,group=input))+geom_line(aes(colour=input),size=1.05)+
    labs(x="Time",y="Percent of HIV +ve population undiagnosed") + 
    geom_vline(xintercept = diag_start,col="midnightblue",size = 0.75) +
    geom_vline(xintercept = art_start,col="midnightblue",size=0.75) + coord_cartesian(ylim = c(0,100))
  
  if(mod_output == FALSE){
  deaths_plot <- ggplot(data = tot_deaths, aes(x = year, y =deaths, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    geom_point(aes(x = year, y = aidsdeaths), colour = "midnightblue", size = 1.5) +
    labs(x = "Time", y= "Number of AIDS deaths")
  
  diagnoses_plots <- ggplot(data = tot_diagnoses, aes(x = year, y = diagnoses, group = input)) +
    geom_line(aes(colour = input), size = 1.05)+
    geom_point(aes(x = year, y = total_cases), size = 1.5, colour = "midnightblue")+
    labs(x = "Time", y = "Diagnoses")
  }else{
    deaths_plot <- ggplot(data = tot_deaths, aes(x = year, y =deaths, group = input)) + 
      geom_line(aes(colour = input), size = 1.05)+
      labs(x = "Time", y= "Number of AIDS deaths")
    
    diagnoses_plots <- ggplot(data = tot_diagnoses, aes(x = year, y = diagnoses, group = input)) +
      geom_line(aes(colour = input), size = 1.05)+
      labs(x = "Time", y = "Diagnoses")
    
  }
  if(mod_output == FALSE){
  diag_rate_plot <- ggplot(data = tot_diag_rate, aes(x = year, y = rate, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    labs(x = "time", y = "rate")
  }else{
    diag_rate_plot <- "nothing to see here"
  }
  art_init_rate_plot <- ggplot(data = art_init_tot, aes(x = time, y = val, group = cat_input)) +
    geom_line(aes(colour = cat, linetype = input), size = 1.01) +
    labs(x = "time", y = "init numbers per class")
  tot_divided_plot <- ggplot(data = tot_divided, aes(x = time, y = deaths, group = class_input)) +
    geom_line(aes(colour = class, linetype = input), size = 1.01) +
    labs(x = "time", y = "deaths from AIDS")
  
  if(mod_output == FALSE){
  combined_plot <- ggpubr::ggarrange(undiag_plot, deaths_plot, diagnoses_plots,
                                     diag_rate_plot, art_init_rate_plot, tot_divided_plot,
                                     ncol = 2, nrow = 3, align = c("v")) 
  }else{
    combined_plot <- ggpubr::ggarrange(undiag_plot, deaths_plot, diagnoses_plots,
                                       art_init_rate_plot, tot_divided_plot,
                                       ncol = 2, nrow = 3, align = c("v")) 
  }
  
  return(list(undiag_df = tot_undiag_data,undiag_plot = undiag_plot, diagnoses_df = tot_diagnoses, diagnoses_plot = diagnoses_plots,
              deaths_df = tot_deaths, deaths_plot = deaths_plot, rate_plot = diag_rate_plot, combined_plot = combined_plot,
              art_inits = art_init_rate_plot, deaths_div = tot_divided_plot))
  
}

testo_lesto <- list(test_run_opt_1, test_run_opt_4)

test_run_plot <- plot_undiagnosed(testo_lesto, mod_output = TRUE)
test_run_plot$combined_plot
test_run_plot$deaths_plot
test_run_plot$deaths_df

matplot(t(art_cd4_counts), type = "l")
test_run_plot$art_inits
art_cd4_counts




