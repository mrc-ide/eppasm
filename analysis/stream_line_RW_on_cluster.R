#####################################################################################
## Strem line with rw ###############################################################
#####################################################################################
setwd("X:")
options(didehpc.username = "jd2117",didehpc.home = "X:/simpleepp",didehpc.cluster = "fi--didemrchnb")

didehpc::didehpc_config(cores = 3,parallel = FALSE)
?didehpc::didehpc_config

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
## !!!!!!!!!!!!!!!!!!!!!!!! Remember to turn on pulse secure at this point !!!!!!!!!!
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
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

#######################################################################################################
## Lets load up and prepare the fp file ###############################################################
#######################################################################################################

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
brazil_fp$t_diagn_start <- 11L   # assume diagnoses starts in 1980


brazil_fp$relinfectART <- 0.3
brazil_fp$tsEpidemicStart <- 1970.5
brazil_fp$likelihood_cd4 <- TRUE
brazil_fp$artinit_use <- FALSE
brazil_fp$aidsdeaths <- TRUE
brazil_fp$time_at_which_get_cd4_counts <- 2001L
brazil_fp$diagn_rate <- array(0.8, c(dim(brazil_fp$cd4_mort), brazil_fp$ss$PROJ_YEARS))
brazil_fp$stages <- 4L

brazil <- list(fp = brazil_fp, csavrd = brazil_csvard_test)

#########################################################################################################
## Now lets get the right Spectrum ART data in there ####################################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/spectrum_ART_numbers",verbose = T)

brazil$fp$art15plus_num <- spectrum_art_numbers

#########################################################################################################
## NOw lets get the corrected ART mortality numbers on there ############################################
#########################################################################################################

load("C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/2015_vals_corrected_ART_mort_brazil", verbose = T)

brazil$fp$art_mort <- test_art_mort

#########################################################################################################
## Now lets get the diagnosis model in there ############################################################
#########################################################################################################

brazil$fp$linear_diagnosis <- "knot_linear"

##########################################################################################################
## Now we can get the art cd4 values in if we want to fit with these, (R MODEL ONLY)!! ###################
##########################################################################################################
load("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/cd4_at_ART_MATRIX", verbose = T)

brazil$fp$tara_art_cd4 <- art_cd4_counts
brazil$fp$use_cd4_art <- TRUE
brazil$fp$time_cd4_stop <- 2015L
brazil$fp$artcd4elig_idx <- rep(1L, 52)

##########################################################################################################
## Lets get our function for plotting the resulting trends ###############################################
##########################################################################################################

plot_undiagnosed <- function(optim_output,diag_start = 1980, art_start = 1996, model_labs){
  tot_undiag_data <- NULL
  input_label <- model_labs
  tot_deaths <- NULL
  tot_diagnoses <- NULL
  tot_diag_rate <- NULL
  art_init_tot <- NULL
  tot_divided <- NULL
  tot_incid <- NULL
  
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
    
    ## INcidence 
    
    incid_df <- cbind.data.frame(list_version$incid15to49, rep(input_label[i], 52), c(1970:2021))
    tot_incid <- rbind.data.frame(tot_incid, incid_df)
    
    
    
  }
  
  names(tot_undiag_data) <- c("undiagnosed","input","year")
  names(tot_deaths) <- c("deaths", "input", "year")
  names(tot_diagnoses) <-c("diagnoses", "input", "year")
  names(tot_diag_rate) <-c("rate","input","year")
  names(art_init_tot) <- c("time","cat","val","input")
  names(tot_divided) <- c("time","class","deaths","input")
  names(tot_incid) <- c("incidence","input","time")
  
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
  
  tot_incid_plot <- ggplot(data = tot_incid, aes(x = time, y = incidence, group = input)) +
    geom_line(aes(colour = input), size = 1.02) +
    labs(x = "time", y = "incidence")
  
  combined_plot <- ggpubr::ggarrange(undiag_plot, deaths_plot, diagnoses_plots,
                                     tot_incid_plot, art_init_rate_plot, tot_divided_plot,
                                     ncol = 2, nrow = 3, align = c("v"), labels = c("Undiagnosed (%)",
                                                                                    "AIDS deaths", "New diagnoses",
                                                                                    "Incidence", "ART initation by CD4 class",
                                                                                    "Deaths by treatment")) 
  
  return(list(undiag_df = tot_undiag_data,undiag_plot = undiag_plot, diagnoses_df = tot_diagnoses, diagnoses_plot = diagnoses_plots,
              deaths_df = tot_deaths, deaths_plot = deaths_plot, rate_plot = diag_rate_plot, combined_plot = combined_plot,
              art_inits = art_init_rate_plot, deaths_div = tot_divided_plot, incid_plot = tot_incid_plot))
  
}

##########################################################################################################
## Now let's test this out for the RW model ##############################################################
##########################################################################################################
brazil$fp$diagnoses_uses <- TRUE
brazil$fp$tARTstart <- 28L

save(brazil, file = "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/updated_brazil_fp_csavrd")
devtools::load_all("C:/Users/josh/Dropbox/hiv_project/eppasm")
devtools::build("C:/Users/josh/Dropbox/hiv_project/eppasm")

brazil_opt_RW <- obj$enqueue(fitmod_csavr(brazil, eppmod="logrw", B0=1e4, optfit=TRUE),
                                         name = "logRW_model")
brazil_opt_RW$status()                                         
brazil_opt_RW$log()
brazilo_opt_log_r <- obj$enqueue(fitmod_csavr(brazil, eppmod="rlogistic", B0=1e4, optfit=TRUE),
                                 name = "optim_knot_linear_mod_3")
brazilo_opt_log_r$status()
brazilo_opt_log_r$log()
brazilo_opt_dooblay_incid <- obj$enqueue(fitmod_csavr(brazil, incid_func = "idbllogistic", B0 = 1e4, optfit = TRUE),
                                         name = "incid double logistic")
brazilo_opt_dooblay_incid$status()
brazilo_opt_dooblay_incid$log()

brazil_r_log_RW <- obj$enqueue(fitmod_csavr(brazil, eppmod = "rlogrw", B0 = 1e4, optfit = TRUE),
                               name = "r log RW")

brazil_r_log_RW$status()
brazil_r_log_RW$log()
opt_rw_res <- brazil_opt_RW$result()
opt_r_log_res <- brazilo_opt_log_r$result()
opt_dooblay_incid <- brazilo_opt_dooblay_incid$result()

opt_res_list <- list(opt_rw_res, opt_r_log_res, opt_dooblay_incid)

output_test <- plot_undiagnosed(opt_res_list, model_labs = c("RW", "rlog","doubleincid"))

output_test$combined_plot

output_test$rate_plot


############################################################################################
## 