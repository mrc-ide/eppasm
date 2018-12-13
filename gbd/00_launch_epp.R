### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm-1/")
date <- substr(gsub("-","",Sys.Date()),3,8)

## Packages
library(data.table)

## Arguments
run.name <- "181126_test"
proj.end <- 2019
n.draws <- 1
cluster.project <- "proj_hiv"

### Paths
dir <- paste0("/ihme/hiv/epp_output/gbd19/", run.name, "/")
dir.create(dir, showWarnings = F)

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
epp.list <- sort(loc.table[epp == 1, ihme_loc_id])
loc.list <- epp.list

# Cache inputs
if(!file.exists(paste0(dir, "populations/"))) {
  prep.job <- paste0("qsub -N eppasm_prep_inputs_", run.name," -P ",cluster.project," -pe multi_slot 5 ",
                      "-e /share/temp/sgeoutput/", user, "/errors ",
                      "-o /share/temp/sgeoutput/", user, "/output ",
                      "/homes/", user, "/HIV/singR_shell.sh ", 
                      code.dir, "R/gbd_prep_inputs.R"," ",run.name," ",proj.end)
  print(prep.job)
  system(prep.job)
}

# Cache prevalence surveys
if(!file.exists(paste0(dir, 'prev_surveys.csv'))){
  prev.job <- paste0("qsub -N eppasm_prev_cache_", run.name," -P ",cluster.project," -pe multi_slot 2 ",
                     "-e /share/temp/sgeoutput/", user, "/errors ",
                     "-o /share/temp/sgeoutput/", user, "/output ",
                     "/homes/", user, "/HIV/singR_shell.sh ", 
                     code.dir, "R/cache_prev_surveys.R"," ",run.name)
  print(prev.job)
  system(prev.job)
}

# Prepare ART proportions
## TODO
# system(paste0("qsub -P ", cluster.project," -N art_prop -hold_jid prev_cache ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/prep_art_props.R ", run.name))

## Launch EPP
for(loc in loc.list) {
    epp.string <- paste0("qsub -P ", cluster.project, " -pe multi_slot 1 ", 
                         "-e /share/temp/sgeoutput/", user, "/errors ",
                         "-o /share/temp/sgeoutput/", user, "/output ",
                         "-N ", loc, "_epp ",
                         "-t 1:", n.draws, " ",
                         "/homes/", user, "/HIV/singR_shell.sh ", 
                         code.dir, "gbd/main.R ",
                         run.name, " ", loc, " ", proj.end)
    print(epp.string)
    system(epp.string)
        
    # Draw compilation
    # draw.string <- paste0("qsub -P ", cluster.project, " -pe multi_slot 2 ", 
    #                       "-e /share/temp/sgeoutput/", user, "/errors ",
    #                       "-o /share/temp/sgeoutput/", user, "/output ",
    #                       "-N ", loc, "_save_draws ",
    #                       "-hold_jid ", loc, "_epp ",
    #                       code.dir, "shell_R.sh ",
    #                       code.dir, "epp-feature-reset/gbd/save_paired_draws.R ",
    #                       loc, " ", run.name, " ", n.draws)
    # print(draw.string)
    # system(draw.string)
}

# ## Check for which draws are missing
# dir <- paste0("/ihme/hiv/epp_output/gbd17/", run.name, "/")
# missing.dt <- rbindlist(lapply(loc.list, function(loc) {
#     files <- list.files(paste0(dir, loc), "incid")
#     numbers <- as.integer(gsub("results_incid", "", gsub(".csv", "", files)))
#     missing <- setdiff(1:n.draws, numbers)
#     if(length(missing) > 0) {
#         dt <- data.table(ihme_loc_id = loc, missing = missing)
#     } else {
#         dt <- data.table()
#     }
# }))
# sum.dt <- missing.dt[, .(num = .N), by = "ihme_loc_id"]
# as.data.frame(sum.dt)
# # loc.list <- sum.dt[num == n.draws, ihme_loc_id]
# 
# ## Fill in missing India and Kenya locations
# system(paste0("qsub -P ", cluster.project," -pe multi_slot 1 -N missing_subnats ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/create_missing_subnats.R ", run.name))
# 
# ## Split India states to Urban Rural
# system(paste0("qsub -P ", cluster.project," -pe multi_slot 1 -N ind_split -hold_jid missing_subnats ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/split_ind_states.R ", run.name))
# 
# ## Plot and combine results
# if(plot) {
#     ind.locs <- c("IND_44538", loc.table[grepl("IND", ihme_loc_id) & level == 5, ihme_loc_id])
#     plot.list <- c(epp.list, ind.locs)
#     for(loc in plot.list) {
#         # Plot individual locations
#         plot.string <- paste0("qsub -P ", cluster.project, " -pe multi_slot 1 ", 
#                               "-e /share/temp/sgeoutput/", user, "/errors ",
#                               "-o /share/temp/sgeoutput/", user, "/output ",
#                               "-N ", loc, "_plot_epp ",
#                               "-hold_jid ", loc, "_save_draws ",
#                               code.dir, "shell_R.sh ",
#                               code.dir, "epp-feature-reset/gbd/plot_epp.R ",
#                               loc, " ", run.name)
#         print(plot.string)
#         system(plot.string)
#     }
# 
#     # Compile all locations
#     combine.holds <- paste(paste0(plot.list, "_plot_epp"), collapse = ",")
#     combine.string <- paste0("qsub -P ", cluster.project, " -pe multi_slot 1 ", 
#                              "-e /share/temp/sgeoutput/", user, "/errors ",
#                              "-o /share/temp/sgeoutput/", user, "/output ",
#                              "-N combine_epp_plots ",
#                              "-hold_jid ", combine.holds, " ",
#                              code.dir, "shell_R.sh ",
#                              code.dir, "epp-feature-reset/gbd/combine_epp_plots.R ",
#                              run.name)
#     print(combine.string)
#     system(combine.string)
# }
# 
# ## Compile comparison data for Spectrum plots
# system(paste0("qsub -P ", cluster.project," -N compile_epp -hold_jid combine_epp_plots ", code.dir, "shell_R.sh ", code.dir, "epp-feature-reset/gbd/compile_epp_data.R ", run.name))

### End
