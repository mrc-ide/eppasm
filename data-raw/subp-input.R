library(tidyverse)
devtools::load_all()

load(file.path(here::here(), "R/sysdata.rda"))
sysdata_obj <- ls()

files <- list.files("~/Documents/Data/Spectrum files/2019 final/subp", "subp$", full.names = TRUE)
subp <- setNames(lapply(files, read_subp_file), sub("^(.*_.*)_.*", "\\1", basename(files)))

subp_gfr2019 <- read.csv("~/Documents/Data/Spectrum files/2019 final/subp_gfr.csv") %>%
  mutate(cc = as.character(cc))

## Keep Tanzania U/R in GFR file
## Add Kavango East GFR; assume same as Kavango
## Tanzania: 'Dodoma' region is misspelled 'Dodmao' --> add a 'Dodmao' entry to dataset

subp_gfr <- bind_rows(
  subp_gfr2019,
  subp_gfr %>% anti_join(subp_gfr2019 %>% select(cc, eppregion)),
  subp_gfr %>% filter(country == "Namibia", eppregion == "Kavango") %>% mutate(eppregion = "Kavango East"),
  subp_gfr %>% filter(country == "Namibia", eppregion == "Kavango") %>% mutate(eppregion = "Kavango West"),
  subp_gfr %>% filter(country == "United Republic of Tanzania", eppregion == "Dodoma") %>% mutate(eppregion = "Dodmoa")
) %>%
  unique 

save(list = sysdata_obj, file = file.path(here::here(), "R/sysdata.rda"))
