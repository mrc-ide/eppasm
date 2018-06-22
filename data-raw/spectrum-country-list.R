
## Country list taken from Spectrum 5.62beta15
## C:\Program Files (x86)\Spectrum5\DP\ModData\CountryListMaster.csv

load("../R/sysdata.rda")

spectrum5_countrylist <- read.csv("CountryListMaster.csv", as.is=TRUE, encoding = "UTF-8")
spectrum5_countrylist$Country[spectrum5_countrylist$Code == 384] <- "Côte d'Ivoire"
devtools::use_data(spectrum5_countrylist, subp, subp_gfr, internal = TRUE, overwrite=TRUE)
