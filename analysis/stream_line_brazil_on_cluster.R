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
                                                                          "C:/Users/josh/Dropbox/hiv_project/epp",
                                                                          "C:/Users/josh/Dropbox/hiv_project/anclik/anclik")))
config <- didehpc::didehpc_config(cores = 2, parallel = FALSE)
obj <- didehpc::queue_didehpc(ctx,config)


######################################## Refresh only !!!!! #######################################
obj$provision(refresh_drat = TRUE)

a <- obj$enqueue(sessionInfo())
a$status()
a$result()
##########################################################################################################################
### NOw lets load up the data to put into the cluster run ################################################################
##########################################################################################################################

library(epp)

library(magrittr)
library(broom)
library(ggplot2)
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/my_created_csvard_for_brazil_WITH_ART",verbose = T)

############################################################################################################
## NOw lets run the thing !!! ##############################################################################
############################################################################################################

brazil_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pjnz)
brazil_fp$artmx_timerr <- rep(1.0, brazil_fp$ss$PROJ_YEARS)
brazil_fp$t_diagn_start <- 11L   # assume diagnoses starts in 1985


brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- TRUE
brazil_fp$artinit_use <- FALSE
brazil_fp$time_at_which_get_cd4_counts <- 2001L

brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$likelihood_cd4 <- TRUE

brazil_fp$stages <- 4L

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

#########################################################################################################
## Now lets get the right ART data in there #############################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/corrected_art_starters_BRAZIL",verbose = T)

brazil$fp$art15plus_num <- art_format

#########################################################################################################
## NOw lets get the corrected ART mortality numbers on there ############################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/2015_vals_corrected_ART_mort_brazil", verbose = T)

brazil$fp$art_mort <- test_art_mort

#########################################################################################################
## Now lets run the model with the raw diagnoses and AIDS deaths and see the output #####################
#########################################################################################################
brazil$fp$aidsdeaths <- TRUE
brazil$fp$diagnoses_uses <- TRUE
brazil$fp$tARTstart <- 28L
brazil$fp$linear_diagnosis <- "gamma"

brazil_fit1_cd4_gamma <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "logistic_single_CD4_diag_deaths_gamma")
brazil_fit1_cd4_gamma$status()
brazil_fit1_cd4_gamma_id <- brazil_fit1_cd4_gamma$id
save(brazil_fit1_cd4_gamma_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/GAMMA_logistic_tara_ART_and_MORT")


brazil_fit2_cd4_gamma <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "double_log_CD4_diag_deaths_gamma")
brazil_fit2_cd4_gamma$status()
brazil_fit2_cd4_gamma_id <- brazil_fit2_cd4_gamma$id
save(brazil_fit2_cd4_gamma_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/GAMMA_double_logistic_tara_ART_and_MORT")


## fit logistic model for transimssion rate (r(t))
brazil_fit3_cd4_gamma <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "rlog_CD4_diag_deaths_gamma")
brazil_fit3_cd4_gamma$status()
brazil_fit3_cd4_gamma_id <- brazil_fit3_cd4_gamma$id
save(brazil_fit3_cd4_gamma_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/GAMMA_r_logistic_tara_ART_and_MORT")


#####################################################################################
#### lets tidy the results and plot them ############################################
#####################################################################################
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")


path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
clusto <- list.files(path_name, full.names = T)
for(i in 1:length(clusto)){
  load(clusto[[i]], verbose = T)
}

fit_1_obj <- obj$task_get(brazil_fit1_cd4_gamma_id)
fit_2_obj <- obj$task_get(brazil_fit2_cd4_gamma_id)
fit_3_obj <- obj$task_get(brazil_fit3_cd4_gamma_id)

brazil_fit1_cd4 <- fit_1_obj$result()
brazil_fit2_cd4 <- fit_2_obj$result()
brazil_fit3_cd4 <- fit_3_obj$result()

brazil_out1_cd4 <- tidy(brazil_fit1_cd4) %>% data.frame(model = "logistic", .)
brazil_out2_cd4 <- tidy(brazil_fit2_cd4) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4 <- tidy(brazil_fit3_cd4) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4)



gamma_diag_and_deaths_updated_art <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  labs(title = "Gamma diagnosis rate, updated ART mortality, fit with CD4 counts of diagnoses and AIDS deaths")

save(gamma_diag_and_deaths_updated_art,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/PLOT_gamma_diag_updated_ART_mort_deaths_and_diagnoses")
save(brazil_out_cd4,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/DATA_gamma_diag_updated_mort_deaths_and_cd4_diagnoses")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
########################################################################################################################
## Now lets queue our model with updated ART mortality numbers, lets run it with the gamma and the linear diagn ########
########################################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

brazil$fp$linear_diagnosis <- "12 param"

brazil_fit1_updated_art_12_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4,
                                                               B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "logistic_single_CD4_diag_deaths_12_linear")
brazil_fit1_updated_art_12_linear$status()
brazil_fit1_updated_art_12_linear_id <- brazil_fit1_updated_art_12_linear$id
save(brazil_fit1_updated_art_12_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/12LINEAR_updated_ART_and_MORT")


brazil_fit2_updated_art_12_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "double_log_CD4_diag_deaths_12_linear")
brazil_fit2_updated_art_12_linear$status()
brazil_fit2_updated_art_12_linear_id <- brazil_fit2_updated_art_12_linear$id
save(brazil_fit2_updated_art_12_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/12LINEAR_double_logistic_updated_ART_and_MORT")


## fit logistic model for transimssion rate (r(t))
brazil_fit3_updated_art_12_linear <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                               name = "rlog_CD4_diag_deaths_12_linear")
brazil_fit3_updated_art_12_linear$status()
brazil_fit3_updated_art_12_linear_id <- brazil_fit3_updated_art_12_linear$id
save(brazil_fit3_updated_art_12_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/12LINEAR_r_logistic_updated_ART_and_MORT")


#### lets tidy the results and plot them #####

path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
clusto <- list.files(path_name, full.names = T)
for(i in 1:length(clusto)){
  load(clusto[[i]], verbose = T)
}

fit_1_obj <- obj$task_get(brazil_fit1_updated_art_12_linear_id)
fit_2_obj <- obj$task_get(brazil_fit2_updated_art_12_linear_id)
fit_3_obj <- obj$task_get(brazil_fit3_updated_art_12_linear_id)

brazil_fit1_cd4 <- fit_1_obj$result()
brazil_fit2_cd4 <- fit_2_obj$result()
brazil_fit3_cd4 <- fit_3_obj$result()

brazil_out1_cd4 <- tidy(brazil_fit1_cd4) %>% data.frame(model = "logistic", .)
brazil_out2_cd4 <- tidy(brazil_fit2_cd4) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4 <- tidy(brazil_fit3_cd4) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4)



no_cd4_fit <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  labs(title = "12 parameter diagnosis rate, fit with CD4 data on diagnoses and AIDS deaths, Spectrum ART mort numbers")

save(no_cd4_fit, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/PLOT_12_param_linear_diagnosis_rate_spectrum_ART_numbers")
save(brazil_out_cd4, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/DATA_12_param_linear_diagnosis_rate_spectrum_ART_numbers")


################################################################################################################
## NOW FOR JEFF'S 6 KNOT POINT DIAGNOSIS MODEL #################################################################
################################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"

brazil_fit1_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e4,
                                                              B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                 name = "logistic_single_CD4_diag_deaths_knot_linear")
brazil_fit1_updated_art_knot_linear$status()
brazil_fit1_updated_art_knot_linear_id <- brazil_fit1_updated_art_knot_linear$id
save(brazil_fit1_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/KNOT_LINEAR_SIMPLE_updated_ART_and_MORT")


brazil_fit2_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic",
                                                                B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                 name = "double_log_CD4_diag_deaths_knot_linear")
brazil_fit2_updated_art_knot_linear$status()
brazil_fit2_updated_art_knot_linear_id <- brazil_fit2_updated_art_knot_linear$id
save(brazil_fit2_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/KNOT_LINEAR_double_log_updated_ART_and_MORT")


## fit logistic model for transimssion rate (r(t))
brazil_fit3_updated_art_knot_linear <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5),
                                                 name = "rlog_CD4_diag_deaths_knot_linear")
brazil_fit3_updated_art_knot_linear$status()
brazil_fit3_updated_art_knot_linear_id <- brazil_fit3_updated_art_knot_linear$id
save(brazil_fit3_updated_art_knot_linear_id,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/KNOT_LINEAR_r_log_updated_ART_and_MORT")


#### lets tidy the results and plot them #####

path_name <- "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_ids/"
clusto <- list.files(path_name, full.names = T)
for(i in 1:length(clusto)){
  load(clusto[[i]], verbose = T)
}

fit_1_obj <- obj$task_get(brazil_fit1_updated_art_knot_linear_id)
fit_2_obj <- obj$task_get(brazil_fit2_updated_art_knot_linear_id)
fit_3_obj <- obj$task_get(brazil_fit3_updated_art_knot_linear_id)

brazil_fit1_cd4 <- fit_1_obj$result()
brazil_fit2_cd4 <- fit_2_obj$result()
brazil_fit3_cd4 <- fit_3_obj$result()

brazil_out1_cd4 <- tidy(brazil_fit1_cd4) %>% data.frame(model = "logistic", .)
brazil_out2_cd4 <- tidy(brazil_fit2_cd4) %>% data.frame(model = "double logistic", .)
brazil_out3_cd4 <- tidy(brazil_fit3_cd4) %>% data.frame(model = "rlogistic", .)
brazil_out_cd4 <- rbind(brazil_out1_cd4, brazil_out2_cd4, brazil_out3_cd4)



knot_linear_diag_plot <- ggplot(subset(brazil_out_cd4, year %in% 1975:2017), 
                     aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line(size=1.05) + geom_ribbon(alpha=0.2) + 
  facet_wrap(~outcome, scales="free") +
  geom_point(aes(y=lik_data), col="darkred", size=0.5) +
  labs(title = "5 Knot diagnosis rate, updated ART mortality rates, fit with AIDS deaths and HIV diagnoses by CD4 count")

save(knot_linear_diag_plot, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/PLOT_knot_linear_diagnosis_rate_updated_ART_mort")
save(brazil_out_cd4, file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/diag_and_deaths_cd4/DATA_knot_linear_diagnosis_rate_updated_ART_mort")

####################################################################################################################
## Working out getting the diagnosis rate and proportion undiagnosed out of these results ##########################
####################################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"
brazil$fp$t_diagn_start <- 11L


brazil_opt1_knot_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),
                                       name = "optim_knot_linear_mod_1")
brazil_opt1_knot_linear$status()
brazil_opt1_knot_linear$log()

brazil_opt2_knot_linear <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_2")
brazil_opt2_knot_linear$status()
brazil_opt2_knot_linear$log()

brazil_opt3_knot_linear <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_3")
brazil_opt3_knot_linear$status()
brazil_opt3_knot_linear$log()

opt_1_res <- brazil_opt1_knot_linear$result()
opt_2_res <- brazil_opt2_knot_linear$result()
opt_3_res <- brazil_opt3_knot_linear$result()

save(opt_2_res,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_knotlinear_mod_2")
save(opt_3_res,
     file = "C:/Users/josh/Dropbox/hiv_project/EPPASM_runs/cluster_runs/cluster_results/optim_knotlinear_mod_3")


opt_list <- list(opt_1_res, opt_2_res, opt_3_res)

plot_undiagnosed <- function(optim_output,diag_start = 1980, art_start = 1996){
  tot_undiag_data <- NULL
  input_label <- paste("input",(1:length(optim_output)),sep = " ")
  tot_deaths <- NULL
  tot_diagnoses <- NULL
  tot_diag_rate <- NULL
  art_init_tot <- NULL
  tot_divided <- NULL
  
  
  for(i in 1:length(optim_output)){
    if(length(optim_output) == 1){
      list_version <- attributes(optim_output[[1]]$mod)
    }else{
      list_version <- attributes(optim_output[[i]]$mod)
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
    
    diag_rate <- optim_output[[i]]$par[(length(optim_output[[i]]$par) - 4) : length(optim_output[[i]]$par)]
    knots <- c(1986, 1996, 2000, 2009, 2015)
    
    diagn_trend <- approx(knots, diag_rate, 1970:2021, rule = 2)$y
    diagn_df <- cbind.data.frame(diagn_trend, rep(input_label[i], 52), c(1970:2021))
    
    tot_diag_rate <- rbind.data.frame(tot_diag_rate, diagn_df)
    
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
  names(tot_diag_rate) <-c("rate","input","year")
  names(art_init_tot) <- c("time","cat","val","input")
  names(tot_divided) <- c("time","class","deaths","input")
  
  art_init_tot$cat_input <- paste(art_init_tot$cat, art_init_tot$input, sep = "")
  tot_divided$class_input <- paste(tot_divided$class, tot_divided$input, sep ="")
  
  
  
  
  tot_deaths <- merge(tot_deaths, optim_output[[1]]$likdat$aidsdeaths, all = TRUE)
  tot_diagnoses <- merge(tot_diagnoses, optim_output[[1]]$likdat$diagnoses, all = TRUE)
  
  undiag_plot <- ggplot(data = tot_undiag_data,aes(x=year,y=undiagnosed,group=input))+geom_line(aes(colour=input),size=1.05)+
    labs(x="Time",y="Percent of HIV +ve population undiagnosed") + 
    geom_vline(xintercept = diag_start,col="midnightblue",size = 0.75) +
    geom_vline(xintercept = art_start,col="midnightblue",size=0.75) + coord_cartesian(ylim = c(0,100))
  
  deaths_plot <- ggplot(data = tot_deaths, aes(x = year, y =deaths, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    geom_point(aes(x = year, y = aidsdeaths), colour = "midnightblue", size = 1.5) +
    labs(x = "Time", y= "Number of AIDS deaths")
  
  diagnoses_plots <- ggplot(data = tot_diagnoses, aes(x = year, y = diagnoses, group = input)) +
    geom_line(aes(colour = input), size = 1.05)+
  geom_point(aes(x = year, y = total_cases), size = 1.5, colour = "midnightblue")+
  labs(x = "Time", y = "Diagnoses")
  
  diag_rate_plot <- ggplot(data = tot_diag_rate, aes(x = year, y = rate, group = input)) + geom_line(aes(colour = input), size = 1.05)+
    labs(x = "time", y = "rate")
  
  art_init_rate_plot <- ggplot(data = art_init_tot, aes(x = time, y = val, group = cat_input)) +
                                 geom_line(aes(colour = cat, linetype = input), size = 1.01) +
                                 labs(x = "time", y = "init numbers per class")
  tot_divided_plot <- ggplot(data = tot_divided, aes(x = time, y = deaths, group = class_input)) +
    geom_line(aes(colour = class, linetype = input), size = 1.01) +
    labs(x = "time", y = "deaths from AIDS")
  
  
  combined_plot <- ggpubr::ggarrange(undiag_plot, deaths_plot, diagnoses_plots,
                                     diag_rate_plot, art_init_rate_plot, tot_divided_plot,
                                     ncol = 2, nrow = 3, align = c("v")) 
  
  return(list(undiag_df = tot_undiag_data,undiag_plot = undiag_plot, diagnoses_df = tot_diagnoses, diagnoses_plot = diagnoses_plots,
              deaths_df = tot_deaths, deaths_plot = deaths_plot, rate_plot = diag_rate_plot, combined_plot = combined_plot,
              art_inits = art_init_rate_plot, deaths_div = tot_divided_plot))
  
}

undiag <- plot_undiagnosed(opt_list)
undiag$combined_plot
ggpubr::annotate_figure(undiag$combined_plot,
                        top = ggpubr::text_grob("Knot_linear diagnosis rates,updated ART mortality, otherwise constant Spectrum values",
                                                color = "red", size = 14))


undiag$diagnoses_plot
undiag$rate_plot

undiag$diagnoses_df

undiag$undiag_plot
undiag$undiag_df
brazil$fp$t_diagn_start
brazil$fp$tARTstart

opt_1_altered_med_cd4 <- brazil_opt1_knot_linear$result()
opt_2_altered_med_cd4 <- brazil_opt2_knot_linear$result()
opt_3_altered_med_cd4 <- brazil_opt3_knot_linear$result()

opt_list_altered_med_cd4 <- list(opt_1_altered_med_cd4, opt_2_altered_med_cd4, opt_3_altered_med_cd4)

undiag_alt <- plot_undiagnosed(opt_list_altered_med_cd4)
undiag_alt$undiag_plot
undiag_alt$diagnoses_plot
undiag_alt$deaths_plot

##############################################################################
## run with the med cat eliminated see what that does ########################
##############################################################################

brazil$fp$artcd4elig_idx <- rep(1L, 52)

brazil_opt1_knot_linear_altered_cd4_index <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),
                                       name = "optim_knot_linear_mod_1")
brazil_opt1_knot_linear_altered_cd4_index$status()
brazil_opt1_knot_linear_altered_cd4_index$log()

brazil_opt2_knot_linear_altered_cd4_index <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_2")
brazil_opt2_knot_linear_altered_cd4_index$status()
brazil_opt2_knot_linear_altered_cd4_index$log()

brazil_opt3_knot_linear_altered_cd4_index <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_3")
brazil_opt3_knot_linear_altered_cd4_index$status()
brazil_opt3_knot_linear_altered_cd4_index$log()

opt_1_res_altered_cd4_index <- brazil_opt1_knot_linear_altered_cd4_index$result()
opt_2_res_altered_cd4_index <- brazil_opt2_knot_linear_altered_cd4_index$result()
opt_3_res_altered_cd4_index <- brazil_opt3_knot_linear_altered_cd4_index$result()

opt_list_altered_cd4_index <- list(opt_1_res_altered_cd4_index,
                                   opt_1_res_altered_cd4_index,
                                   opt_3_res_altered_cd4_index)

undiag__altered_cd4_index <- plot_undiagnosed(opt_list_altered_cd4_index)


undiag_eliminated$undiag_plot
undiag_eliminated$diagnoses_plot
undiag_eliminated$deaths_plot
undiag_eliminated$rate_plot

undiag__altered_cd4_index$combined_plot
ggpubr::annotate_figure(undiag__altered_cd4_index$combined_plot,
                        top = ggpubr::text_grob("Knot_linear diagnosis rates,updated ART mortality, ART eligibility index set to 1",
                                                color = "red", size = 14))

#############################################################################################################
## Now lets try altering the median_cd4init_ all to 0, as in 2009 this has some value #######################
#############################################################################################################
brazil$fp$linear_diagnosis <- "knot_linear"
brazil$fp$t_diagn_start <- 11L

which(brazil$fp$median_cd4init != 0)

brazil$fp$median_cd4init[which(brazil$fp$median_cd4init != 0)] <- 0L
brazil$fp$med_cd4init_cat < - rep(0L, 52)
brazil$fp$artcd4elig_idx <- brazil_fp$artcd4elig_idx

brazil_opt1_median_cd4_init_to_0 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),
                                       name = "optim_knot_linear_mod_1")
brazil_opt1_median_cd4_init_to_0$status()
brazil_opt1_median_cd4_init_to_0$log()

brazil_opt2_median_cd4_init_to_0 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_2")
brazil_opt2_median_cd4_init_to_0$status()
brazil_opt2_median_cd4_init_to_0$log()

brazil_opt3_median_cd4_init_to_0 <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=5e4, optfit=TRUE),
                                       name = "optim_knot_linear_mod_3")
brazil_opt3_median_cd4_init_to_0$status()
brazil_opt3_median_cd4_init_to_0$log()

opt_1_res_eliminated_median <- brazil_opt1_median_cd4_init_to_0$result()
opt_2_res_eliminated_median <- brazil_opt2_median_cd4_init_to_0$result()
opt_3_res_eliminated_median <- brazil_opt3_median_cd4_init_to_0$result()

opt_list_eliminated_median <- list(opt_1_res_eliminated_median, opt_1_res_eliminated_median,
                                    opt_3_res_eliminated_median)

undiag_eliminated_median <- plot_undiagnosed(opt_list_eliminated_median)
undiag_eliminated_median$art_inits
undiag_eliminated_median$combined_plot

ggpubr::annotate_figure(undiag_eliminated_median$combined_plot,
                        top = ggpubr::text_grob("Knot_linear diagnosis rates,updated ART mortality, Median CD4 set to 0 and Med CD4 cat set to 0",
                                                color = "red", size = 14))


##########################################################################################
## Setting all CD4 catergories to baseline ###############################################
##########################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"
brazil$fp$t_diagn_start <- 11L

which(brazil$fp$median_cd4init != 0)

brazil$fp$median_cd4init[which(brazil$fp$median_cd4init != 0)] <- 0L
brazil$fp$med_cd4init_cat <- rep(0L, 52)
brazil$fp$artcd4elig_idx <- rep(1L, 52)
brazil$fp$med_cd4init_input <- rep(0L, 52)


brazil_opt1_all_to_0 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),
                                                name = "optim_knot_linear_mod_1")
brazil_opt1_all_to_0$status()
brazil_opt1_all_to_0$log()

brazil_opt2_all_to_0 <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                                name = "optim_knot_linear_mod_2")
brazil_opt2_all_to_0$status()
brazil_opt2_all_to_0$log()

brazil_opt3_all_to_0 <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                                name = "optim_knot_linear_mod_3")
brazil_opt3_all_to_0$status()
brazil_opt3_all_to_0$log()

opt_1_res_all_to_0 <- brazil_opt1_all_to_0$result()
opt_2_res_all_to_0 <- brazil_opt2_all_to_0$result()
opt_3_res_all_to_0 <- brazil_opt3_all_to_0$result()

opt_list_all_to_0 <- list(opt_1_res_all_to_0, opt_1_res_all_to_0,
                                   opt_3_res_all_to_0)

undiag_all_to_0 <- plot_undiagnosed(opt_list_all_to_0)
undiag_all_to_0$art_inits
undiag_all_to_0$combined_plot

ggpubr::annotate_figure(undiag_all_to_0$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates,updated ART mortality, all CD4 catergories set to 1 or 0",
                                                color = "red", size = 14))



####################################################################################################################
## Artificially setting the eligibility criteria to only those who have aids to see affect on hiv deaths ###########
####################################################################################################################

brazil$fp$median_cd4init[which(brazil$fp$median_cd4init != 0)] <- 100L 
brazil$fp$med_cd4init_cat <- rep(5L, 52)
brazil$fp$med_cd4init_input <- rep(0L, 52)
brazil$fp$artcd4elig_idx <- rep(5L, 52)
brazil$fp$med_cd4init_input <- brazil_fp$med_cd4init_input


brazil_opt1_elig_only_AIDS <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=5e4, optfit=TRUE),
                                    name = "optim_knot_linear_mod_1")
brazil_opt1_elig_only_AIDS$status()
brazil_opt1_elig_only_AIDS$log()

brazil_opt2_elig_only_AIDS <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                    name = "optim_knot_linear_mod_2")
brazil_opt2_elig_only_AIDS$status()
brazil_opt2_elig_only_AIDS$log()

brazil_opt3_elig_only_AIDS <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                    name = "optim_knot_linear_mod_3")
brazil_opt3_elig_only_AIDS$status()
brazil_opt3_elig_only_AIDS$log()

opt_1_elig_only_AIDS <- brazil_opt1_elig_only_AIDS$result()
opt_2_elig_only_AIDS <- brazil_opt2_elig_only_AIDS$result()
opt_3_elig_only_AIDS <- brazil_opt3_elig_only_AIDS$result()

opt_list_elig_only_AIDS <- list(opt_1_elig_only_AIDS, opt_1_elig_only_AIDS,
                          opt_3_elig_only_AIDS)

undiag_elig_only_AIDS <- plot_undiagnosed(opt_list_elig_only_AIDS)
undiag_elig_only_AIDS$art_inits
undiag_elig_only_AIDS$combined_plot

ggpubr::annotate_figure(undiag_elig_only_AIDS$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates,updated ART mortality, only AIDS patients eligible for ART",
                                                color = "red", size = 14))














##### all init cats sets to 0

opt_1_res_eliminated_all <- brazil_opt1_all_0$result()
opt_2_res_eliminated_all <- brazil_opt2_all_0$result()
opt_3_res_eliminated_all <- brazil_opt3_all_0$result()

opt_list_eliminated_all <- list(opt_1_res_eliminated_all, opt_1_res_eliminated_all,
                                   opt_3_res_eliminated_all)

undiag_eliminated_all <- plot_undiagnosed(opt_list_eliminated_all)
undiag_eliminated_all$combined_plot
undiag_eliminated_all$undiag_df

#######################################
## fitting with art inits as well #####
#######################################

opt_1_res_eliminated_median_art_init <- brazil_opt1_all_0_art_init_fit$result()
opt_2_res_eliminated_median_art_init <- brazil_opt2_all_0_art_init_fit$result()
opt_3_res_eliminated_median_art_init <- brazil_opt3_all_0_art_init_fit$result()

opt_list_eliminated_median_art_init <- list(opt_1_res_eliminated_median_art_init, opt_1_res_eliminated_median_art_init,
                                   opt_3_res_eliminated_median_art_init)

undiag_eliminated_median_art_init <- plot_undiagnosed(opt_list_eliminated_median_art_init)
undiag_eliminated_median_art_init$combined_plot


undiag_eliminated_median$art_inits
undiag_eliminated_all$art_inits
undiag_eliminated_median_art_init$art_inits
undiag_eliminated_median_art_init$combined_plot
test_simmod <- simmod.specfp(brazil$fp, VERSION = "R")


#####################################################################
## Using knot_linear with original artcd4 elig index ################
#####################################################################
brazil_opt1_orig_cd4_elig <- obj$enqueue(fitmod_csavr(brazil, incid_func = "ilogistic", B0=1e3, optfit=TRUE),
                                              name = "optim_knot_linear_mod_1")
brazil_opt1_orig_cd4_elig$status()
brazil_opt1_orig_cd4_elig$log()

brazil_opt2_orig_cd4_elig <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                              name = "optim_knot_linear_mod_2")
brazil_opt2_orig_cd4_elig$status()
brazil_opt2_orig_cd4_elig$log()

brazil_opt3_orig_cd4_elig <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                              name = "optim_knot_linear_mod_3")
brazil_opt3_orig_cd4_elig$status()
brazil_opt3__orig_cd4_elig$log()

opt_1_res_orig_cd4_elig <- brazil_opt1_orig_cd4_elig$result()
opt_2_res_orig_cd4_elig <- brazil_opt2_orig_cd4_elig$result()
opt_3_res_orig_cd4_elig <- brazil_opt3_orig_cd4_elig$result()

opt_list_orig_cd4_elig <- list(opt_1_res_orig_cd4_elig, opt_1_res_orig_cd4_elig,
                                   opt_3_res_orig_cd4_elig)

undiag_orig_cd4_elig <- plot_undiagnosed(opt_list_orig_cd4_elig)
undiag_orig_cd4_elig$art_inits
undiag_orig_cd4_elig$combined_plot

save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_everything_spectrum_standard")
brazil$fp$artcd4elig_idx <- rep(1L, 52)
save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_artcd4_elig_index_changed_to_1")

#########################################################################################################
## Using Jeff's Spectrum ART numbers ####################################################################
#########################################################################################################

jeff_numbers <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/adult-art-by-sex_Spectrum2018.csv")
jeff_numbers$year <- seq(1997, by = 1, length.out = nrow(jeff_numbers))

brazil$fp$art15plus_num[1,28:52] <- jeff_numbers$ï..Male
brazil$fp$art15plus_num[2,28:52] <- jeff_numbers$Female

spectrum_art_numbers <- brazil$fp$art15plus_num

save(spectrum_art_numbers,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/spectrum_ART_numbers")

#### Run with even elgibility first ####

brazil_opt2_even_elig <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                         name = "optim_knot_linear_mod_2")
brazil_opt2_even_elig$status()
brazil_opt2_even_elig$log()

brazil_opt3_even_elig <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                         name = "optim_knot_linear_mod_3")
brazil_opt3_even_elig$status()
brazil_opt3_even_elig$log()

opt_2_res_even_elig <- brazil_opt2_even_elig$result()
opt_3_res_even_elig <- brazil_opt3_even_elig$result()

opt_list_even_elig <- list(opt_2_res_even_elig,
                               opt_3_res_even_elig)

undiag_even_elig <- plot_undiagnosed(opt_list_even_elig)
undiag_even_elig$art_inits
undiag_even_elig$combined_plot

ggpubr::annotate_figure(undiag_even_elig$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates, all CD4 classes eligible  Spectrum ART numbers",
                                                color = "red", size = 14))

save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_everything_spectrum_standard")
brazil$fp$artcd4elig_idx <- rep(1L, 52)
save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_artcd4_elig_index_changed_to_1")

######################
## orig cd4 elig #####
######################

brazil$fp$artcd4elig_idx <- brazil_fp$artcd4elig_idx 

brazil_opt2_orig_elig <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0=1e4, optfit=TRUE),
                                     name = "optim_knot_linear_mod_2")
brazil_opt2_orig_elig$status()
brazil_opt2_orig_elig$log()

brazil_opt3_orig_elig <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                     name = "optim_knot_linear_mod_3")
brazil_opt3_orig_elig$status()
brazil_opt3_orig_elig$log()

opt_2_res_orig_elig <- brazil_opt2_orig_elig$result()
opt_3_res_orig_elig <- brazil_opt3_orig_elig$result()

opt_list_orig_elig <- list(opt_2_res_orig_elig,
                           opt_3_res_orig_elig)

undiag__orig_elig <- plot_undiagnosed(opt_list_orig_elig)
undiag_even_elig$art_inits
undiag__orig_elig$combined_plot

ggpubr::annotate_figure(undiag__orig_elig$combined_plot,
                        top = ggpubr::text_grob("Knot linear diagnosis rates, original spectrum classes, Spectrum ART numbers",
                                                color = "red", size = 14))



save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_everything_spectrum_standard")
brazil$fp$artcd4elig_idx <- rep(1L, 52)
save(brazil,
     file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/brazil_obj_artcd4_elig_index_changed_to_1")



