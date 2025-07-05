
## Country list taken from Spectrum 5.62beta15
## C:\Program Files (x86)\Spectrum5\DP\ModData\CountryListMaster.csv

## 15 June 2024:
## - Manually updated with Dominica and St. Kitts & Nevis

## 5 July 2025:
## - Manually added Tuvalu to country list


load("../R/sysdata.rda")

spectrum5_countrylist <- utils::read.csv("CountryListMaster_2025-07-05.csv",
                                         as.is=TRUE, encoding = "UTF-8")
spectrum5_countrylist$Country[spectrum5_countrylist$Code == 384] <- "Côte d'Ivoire"

iso <- utils::read.csv("iso.csv", na.strings = "", as.is=TRUE)
names(iso) <- sub("alpha.2", "iso2", names(iso))
names(iso) <- sub("alpha.3", "iso3", names(iso))


spectrum5_countrylist <- merge(spectrum5_countrylist, iso[c("country.code", "iso2", "iso3")],
                               by.x="Code", by.y="country.code", all.x=TRUE)

usethis::use_data(spectrum5_countrylist, subp, subp_gfr, internal = TRUE, overwrite=TRUE)
