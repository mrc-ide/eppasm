### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))

## Packages
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
	run.name <- args[1]
	proj.end <- args[2]
} else {
	run.name <- "181126_test"
	proj.end <- 2019
}

out.dir <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/")
dir.create(out.dir, showWarnings = F)

## Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
source(paste0(root, "/temp/central_comp/libraries/current/r/get_population.R"))
source('/home/j/temp/central_comp/libraries/2019_gbd_env/r/get_covariate_estimates.R')

## Locations
loc.table <- get_locations(hiv_metadata = TRUE)
write.csv(loc.table, paste0(out.dir, 'location_table.csv'), row.names = F)
epp.locs <- loc.table[epp == 1, location_id]
parent.locs <- loc.table[(grepl('IND', ihme_loc_id) & level < 5) | (grepl('KEN', ihme_loc_id) & level < 5) | ihme_loc_id %in% c('NGA', 'ETH', 'ZAF'), location_id]

## Population
pop <- get_population(age_group_id = c(28, 49:128), location_id = epp.locs, year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2, single_year_age = T)
dir.create(paste0(out.dir, '/population_single_age'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id) {
  out.pop <- copy(pop[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(out.pop, paste0(out.dir, '/population_single_age/', c.iso, ".csv"), row.names = F)
}))

pop <- get_population(age_group_id = c(8:20), location_id = c(epp.locs, parent.locs), year_id = 1970:2019, gbd_round_id = 6, sex_id = 1:2)
dir.create(paste0(out.dir, '/population'), showWarnings = F)
invisible(lapply(c(epp.locs, parent.locs), function(c.location_id) {
  out.pop <- copy(pop[location_id == c.location_id])
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(out.pop, paste0(out.dir, '/population/', c.iso, ".csv"), row.names = F)
}))

pop.15to49 <- get_population(age_group_id = 8:15, location_id = epp.locs, year_id = 1950:2019, gbd_round_id = 6, sex_id = 1:2)
dir.create(paste0(out.dir, '/population_15to49'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id) {
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(pop.15to49[location_id == c.location_id], paste0(out.dir, '/population_15to49/', c.iso, ".csv"), row.names = F)
}))

## Migration
## TODO: migration filepath needs to be updated for GBD 2019
mig <- fread(paste0('/ihme/fertilitypop/gbd_2017/population/modeling/popReconstruct/v96/best/net_migrants.csv'))
setnames(mig, c('year_id', 'sex_id'), c('year', 'sex'))
mig[age > 80, age := 80]
mig <- mig[year >= 1970, .(value = sum(value)), by = c('age', 'sex', 'year', 'ihme_loc_id')]
dir.create(paste0(out.dir, '/migration'), showWarnings = F)
invisible(lapply(epp.locs, function(c.location_id){
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  mig.loc <- mig[ihme_loc_id == c.iso]
  if(nrow(mig.loc) == 0){
    if(loc.table[ihme_loc_id == c.iso, level] > 3){
      ## Crude population split of migration for subnationals -- this is temporary until pop model outputs subnationals
      parent.iso <- substr(c.iso, 1, 3)
      parent.pop <- fread(paste0(out.dir, '/population/', parent.iso, '.csv'))
      child.pop <- fread(paste0(out.dir, '/population/', c.iso, '.csv'))
      setnames(parent.pop, "population", "parent")
      setnames(child.pop, "population", "child")
      merged.pop <- merge(parent.pop, child.pop, by = c("age_group_id", "year_id", "sex_id"))
      collapsed.pop <- merged.pop[, lapply(.SD, sum), by = .(year_id, sex_id), .SDcols = c("parent", "child")]
      collapsed.pop[, prop := child / parent]
      setnames(collapsed.pop, c('year_id', 'sex_id'), c('year', 'sex'))
      mig.loc <- mig[ihme_loc_id == parent.iso]
      mig.loc <- merge(mig.loc, collapsed.pop, by = c('year', 'sex'))
      mig.loc[, value := value * prop]
    } else{
      mig.loc <- data.table(expand.grid(age = 0:80, sex = 1:2, year = 1970:proj.end, ihme_loc_id = c.iso, value = 0))
    }
  }
  mig.loc[, c('parent', 'child', 'prop', 'ihme_loc_id') := NULL]
  write.csv(mig.loc, paste0(out.dir, '/migration/', c.iso, '.csv'), row.names = F)
}))


## ASFR
asfr <- get_covariate_estimates(covariate_id = 13, location_id = epp.locs)
asfr <- asfr[age_group_id %in% c(8:14) & sex_id == 2, list(year_id, age_group_id, mean_value, location_id)]
asfr[, age := (age_group_id - 5) * 5]
setnames(asfr, c('mean_value', 'year_id'), c('value', 'year'))
dir.create(paste0(out.dir, '/ASFR'))
invisible(lapply(epp.locs, function(c.location_id){
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(asfr[location_id == c.location_id, list(value, age, year)], paste0(out.dir, '/ASFR/', c.iso, '.csv'), row.names = F)
}))

## TFR
tfr <- get_covariate_estimates(covariate_id = 149, location_id = epp.locs)
tfr <- tfr[,list(location_id, year_id, mean_value)]
setnames(tfr, c('mean_value', 'year_id'), c('value', 'year'))
dir.create(paste0(out.dir, '/TFR'))
invisible(lapply(epp.locs, function(c.location_id){
  c.iso <- loc.table[location_id == c.location_id, ihme_loc_id]
  write.csv(tfr[location_id == c.location_id, list(value, year)], paste0(out.dir, '/TFR/', c.iso, '.csv'), row.names = F)
}))
