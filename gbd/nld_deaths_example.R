deaths <- fread('/snfs1/WORK/04_epi/01_database/02_data/hiv/spectrum/summary/180702_numbat_combined/locations/NLD_spectrum_prep.csv')
deaths <- deaths[sex_id %in% 1:2 & age_group_id >= 5 & !age_group_id %in% c(20, 22, 24) & measure == 'Deaths' & metric == 'Count']
deaths <- deaths[,.(mean = sum(mean)), by = c('year_id', 'sex_id')]

head(deaths)

deaths <- reshape2::acast(deaths, sex_id ~ year_id)


deaths_dt <- deaths[ , match(rep(as.character(1971:2017), each = 10), colnames(deaths))] / 10

nl_fp <- prepare_directincid("/homes/tahvif/Netherlands2017.PJNZ")
fp <- nl_fp
mod <- simmod(fp)

fp$deaths_dt <- deaths_dt

match(colnames(deaths))

debugonce(simmod.specfp)

modR <- simmod(fp, "R")

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
