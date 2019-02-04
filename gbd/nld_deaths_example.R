nl_fp <- prepare_directincid("/homes/tahvif/Netherlands2017.PJNZ")
fp <- nl_fp
mod <- simmod(fp)


deaths <- fread('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/180702_numbat_combined/locations/NLD_spectrum_prep.csv')
deaths <- deaths[sex_id %in% 1:2 & age_group_id >= 5 & !age_group_id %in% c(20, 22, 24) & measure == 'Deaths' & metric == 'Count']
deaths <- deaths[,.(mean = sum(mean)), by = c('year_id', 'sex_id')]

head(deaths)
## Placeholder for now
while(length(unique(deaths$year_id)) != fp$SIM_YEARS){
  deaths.ext <- deaths[year_id == max(deaths$year_id)]
  deaths.ext[, year_id := year_id + 1]
  deaths <- rbind(deaths, deaths.ext, use.names = T)
}
deaths <- reshape2::acast(deaths, sex_id ~ year_id)


deaths_dt <- deaths[ , match(rep(as.character(1970:2021), each = 10), colnames(deaths))] / 10

fp$deaths_dt <- deaths_dt
setwd('/homes/tahvif/eppasm')
devtools::load_all()

## Likelihood for VR data
## 1. Calculate expected number of deaths predicted by model
## 2. Calculate log likelihood using Poisson

expected_deaths <- colSums(sweep(attr(mod, "artpop"), 1:4, fp$art_mort, "*"),,3) + 
  colSums(sweep(attr(mod, "hivpop"), 1:3, fp$cd4_mort, "*"),,2)

sum(ldpois(deaths[ , 12:48], expected_deaths[ , 12:48]))

fp$deaths_dt[,ts]
expected_deaths <- DT * cd4_mort_ts * hivpop[,,,i]

delta_ts <- fp$deaths_dt[,ts] / colSums(expected_deaths,,2)

ldpois <- function(x, lambda){
  x * log(lambda) - lambda - lgamma(x+1)
}

grad <- grad - cd4_mort_ts * hivpop[,,,i]              # HIV mortality, untreated


## Tahvi
eppd <- list(vr = deaths, 
             country = 'Netherlands',
             region = 'NLD',
             projset_id = 0)

dt <- list()
fp$ss$time_epi_start <- 1971
fp$incidinput <- NULL
fp$incidpopage <- NULL
attr(dt, 'specfp') <- fp
attr(dt, 'eppd') <- eppd
fit <- fitmod(dt, eppmod = 'rhybrid', rw_start = 2010, B0=1e2, B=1e2, VERSION = 'R', number_k = 10)
gbd_sim_mod(fit, VERSION = 'R')
rand.draw <- round(runif(1, min = 1, max = 3000))
output.dt <- get_gbd_outputs(result[[rand.draw]], attr(dt, 'specfp'))

