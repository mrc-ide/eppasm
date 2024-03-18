###############################################################################
####                                                                       ####
####  This script assembles subnational general fertility rates (GFR) for  ####
####  women aged 15-44 for EPP subnational regions.                        ####
####                                                                       ####
###############################################################################



##  Get DHS GFR estimates from StatCompiler via the DHS API
devtools::load_all("~/Documents/Code/R/hhsurveydata")  # www.github.com/mrc-ide/hhsurveydata
indicators <- jsonlite::fromJSON("https://api.dhsprogram.com/rest/dhs/indicators?apiKey=ICLSPH-527168")$Data

## GFR from most recent survey in each country, with all stratifications
gfr <- fetch_dhsapi("https://api.dhsprogram.com/rest/dhs/data/FE_FRTR_W_GFR,all?selectSurveys=byIndicator&indicatorIds=FE_FRTR_W_GFR&apiKey=ICLSPH-527168")


## EPP-ASM regions
load("~/Documents/Code/R/eppasm/R/sysdata.rda")

eppasm_subp <- mapply(data.frame, file=names(subp), eppregion=lapply(subp, names), SIMPLIFY=FALSE)
eppasm_subp <- do.call(rbind, eppasm_subp)
eppasm_subp$cc <- sub("(.*)_(.*)", "\\1", eppasm_subp$file)
eppasm_subp$country <- sub("(.*)_(.*)", "\\2", eppasm_subp$file)
rownames(eppasm_subp) <- NULL

countries <- sub(".*_(.*)", "\\1", names(subp))
setdiff(countries, gfr$CountryName)

gfr$gfr <- gfr$Value / 1000
gfr$survyear <- gfr$SurveyYear
gfr$survtype <- gfr$SurveyType
gfr$country <- gfr$CountryName
gfr$country <- sub("Tanzania", "United Republic of Tanzania", gfr$country)
gfr$country <- sub("Congo Democratic Republic", "Democratic Republic of the Congo", gfr$country)

gfr$eppregion <- gfr$CharacteristicLabel
gfr$eppregion <- sub("Total", "N", gfr$eppregion)
gfr$eppregion <- sub("Urban", "U", gfr$eppregion)
gfr$eppregion <- sub("Rural", "R", gfr$eppregion)
gfr$eppregion <- sub("Northern", "Northern Region", gfr$eppregion)
gfr$eppregion <- sub("Central", "Central Region", gfr$eppregion)
gfr$eppregion <- sub("Southern", "Southern Region", gfr$eppregion)

gfr$eppregion[gfr$country == "Swaziland" & gfr$CharacteristicCategory == "Region"] <-
  paste(gfr$eppregion[gfr$country == "Swaziland" & gfr$CharacteristicCategory == "Region"], "Region")
gfr$eppregion <- sub("\\.\\.", "", gfr$eppregion)


subp_gfr <- merge(eppasm_subp, gfr[c("country", "eppregion", "survyear", "survtype", "gfr")], all.x=TRUE)

subp_gfr <- subp_gfr[c("file", "cc", "country", "eppregion", "survyear", "survtype", "gfr")]
subp_gfr <- subp_gfr[with(subp_gfr, order(file, cc, country, eppregion)),]



#############################################################
####  GFR estimates not available from DHS StatCompiler  ####
#############################################################

## Cape Verde: DHS 2005 report, Table 4.1 (https://dhsprogram.com/pubs/pdf/FR203/FR203.pdf)
subset(subp_gfr, country == "Cabo Verde")
subp_gfr[subp_gfr$country == "Cabo Verde" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2005, "DHS", 0.098)
subp_gfr[subp_gfr$country == "Cabo Verde" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2005, "DHS", 0.093)
subp_gfr[subp_gfr$country == "Cabo Verde" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2005, "DHS", 0.106)

unique(subset(subp_gfr, is.na(gfr))$country)


## Mauritania: Special DHS 2003-04 report, Table 3.1 (https://dhsprogram.com/pubs/pdf/FR150/FR150.pdf)
subset(subp_gfr, country == "Mauritania")
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2004, "DHS", 0.146)
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2004, "DHS", 0.131)
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2004, "DHS", 0.161)


## Equatorial Guinea: DHS 2011 report, Table 4.1 (https://dhsprogram.com/pubs/pdf/FR271/FR271.pdf)
subset(subp_gfr, country == "Equatorial Guinea")
subp_gfr[subp_gfr$country == "Equatorial Guinea" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2011, "DHS", 0.181)
subp_gfr[subp_gfr$country == "Equatorial Guinea" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2011, "DHS", 0.161)
subp_gfr[subp_gfr$country == "Equatorial Guinea" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2011, "DHS", 0.201)


## South Africa: DHS 2015 Key Findings report, Table 3 (https://dhsprogram.com/pubs/pdf/PR84/PR84.pdf)
subset(subp_gfr, country == "South Africa")
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2015, "DHS", 0.094)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2015, "DHS", 0.087)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2015, "DHS", 0.109)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "EC", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.089)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "FS", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.092)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "GP", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.077)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "KZN", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.090)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "LP", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.107)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "MP", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.092)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "NC", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.094)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "NW", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.096)
subp_gfr[subp_gfr$country == "South Africa" & subp_gfr$eppregion == "WC", c("survyear", "survtype", "gfr")] <- c(2015, "Thembisa", 0.081)

## South Sudan: MICS 2010

wm <- foreign::read.spss("~/Documents/Data/MICS/South Sudan_MICS4_Datasets/South Sudan MICS 2010 SPSS Datasets/wm.sav")
wm.variable.labels <- attr(wm, "variable.labels")
wm.label.tabel <- attr(wm, "label.table")
wm <- data.frame(wm)

wm$start <- pmax(wm$WDOI-36, wm$WDOB+15*12)
wm$end <- pmin(wm$WDOI, wm$WDOB+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$wmweight
wm <- subset(wm, pys_wt > 0)


bh <- foreign::read.spss("~/Documents/Data/MICS/South Sudan_MICS4_Datasets/South Sudan MICS 2010 SPSS Datasets/bh.sav")
bh.variable.labels <- attr(bh, "variable.labels")
bh.label.tabel <- attr(bh, "label.table")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, BH4C > start & BH4C <= end)
bh$births_wt <- bh$wmweight

gfr <- merge(stats::aggregate(births_wt ~ HH6, bh, sum),
             stats::aggregate(pys_wt ~ HH6, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)

subset(subp_gfr, country == "South Sudan")
subp_gfr[subp_gfr$country == "South Sudan" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2010, "MICS", 0.237)
subp_gfr[subp_gfr$country == "South Sudan" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2010, "MICS", 0.237)
subp_gfr[subp_gfr$country == "South Sudan" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2010, "MICS", 0.237)



## Guinea-Bissau: MICS 2014

wm <- foreign::read.spss("~/Documents/Data/MICS/Guinea Bissau_MICS5_Datasets/Guinea Bissau MICS 2014 SPSS Datasets/wm.sav")
wm.variable.labels <- attr(wm, "variable.labels")
wm.label.tabel <- attr(wm, "label.table")
wm <- data.frame(wm)

wm$start <- pmax(wm$WDOI-36, wm$WDOB+15*12)
wm$end <- pmin(wm$WDOI, wm$WDOB+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$wmweight
wm <- subset(wm, pys_wt > 0)


bh <- foreign::read.spss("~/Documents/Data/MICS/Guinea Bissau_MICS5_Datasets/Guinea Bissau MICS 2014 SPSS Datasets/bh.sav")
bh.variable.labels <- attr(bh, "variable.labels")
bh.label.tabel <- attr(bh, "label.table")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, BH4C > start & BH4C <= end)
bh$births_wt <- bh$wmweight

gfr <- merge(stats::aggregate(births_wt ~ HH6, bh, sum),
             stats::aggregate(pys_wt ~ HH6, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)

subset(subp_gfr, country == "Guinea-Bissau")
subp_gfr[subp_gfr$country == "Guinea-Bissau" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.159)
subp_gfr[subp_gfr$country == "Guinea-Bissau" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.114)
subp_gfr[subp_gfr$country == "Guinea-Bissau" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.204)


## Mauritania: MICS 2011

wm <- foreign::read.spss("~/Documents/Data/MICS/Mauritania_MICS4_Datasets/Mauritania MICS 2011 SPSS Datasets/wm.sav")
wm <- data.frame(wm)

wm$start <- pmax(wm$WDOI-36, wm$WDOB+15*12)
wm$end <- pmin(wm$WDOI, wm$WDOB+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$WMWEIGHT
wm <- subset(wm, pys_wt > 0)


bh <- foreign::read.spss("~/Documents/Data/MICS/Mauritania_MICS4_Datasets/Mauritania MICS 2011 SPSS Datasets/bh.sav")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, BH4C > start & BH4C <= end)
bh$births_wt <- bh$WMWEIGHT

gfr <- merge(stats::aggregate(births_wt ~ HH6, bh, sum),
             stats::aggregate(pys_wt ~ HH6, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)
sum(gfr$births_wt) / sum(gfr$pys_wt)

subset(subp_gfr, country == "Mauritania")
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2011, "MICS", 0.159)
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2011, "MICS", 0.139)
subp_gfr[subp_gfr$country == "Mauritania" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2011, "MICS", 0.177)


## Swaziland: MICS 2014

wm <- foreign::read.spss("~/Documents/Data/MICS/Swaziland_MICS5_Datasets/Swaziland MICS 2014 SPSS Datasets/wm.sav")
wm <- data.frame(wm)

wm$start <- pmax(wm$WDOI-36, wm$WDOB+15*12)
wm$end <- pmin(wm$WDOI, wm$WDOB+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$wmweight
wm <- subset(wm, pys_wt > 0)


bh <- foreign::read.spss("~/Documents/Data/MICS/Swaziland_MICS5_Datasets/Swaziland MICS 2014 SPSS Datasets/bh.sav")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, BH4C > start & BH4C <= end)
bh$births_wt <- bh$wmweight

gfr <- merge(stats::aggregate(births_wt ~ HH6, bh, sum),
             stats::aggregate(pys_wt ~ HH6, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)
gfr
sum(gfr$births_wt) / sum(gfr$pys_wt)

subset(subp_gfr, country == "Swaziland")
subp_gfr[subp_gfr$country == "Swaziland" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.117)
subp_gfr[subp_gfr$country == "Swaziland" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.106)
subp_gfr[subp_gfr$country == "Swaziland" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.122)


## Somalia: MICS 2006

wm <- foreign::read.spss("~/Documents/Data/MICS/Somalia 2006 MICS_Datasets/Somalia MICS 2006 SPSS Datasets/wm.sav")
wm <- data.frame(wm)

wm$start <- pmax(wm$cmcdoiw-36, wm$wdob+15*12)
wm$end <- pmin(wm$cmcdoiw, wm$wdob+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$wmweight
wm <- subset(wm, pys_wt > 0)

wm$HH6 <- factor(wm$HH6, c("Urban", "Non Urban"), c("Urban", "Rural"))
wm$HH7 <- factor(wm$HH7, c("North West", "North East", "Central South"), c("Somalialand", "Puntland", "Central South"))


bh <- foreign::read.spss("~/Documents/Data/MICS/Somalia 2006 MICS_Datasets/Somalia MICS 2006 SPSS Datasets/bh.sav")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, IMPBH4C > start & IMPBH4C <= end)
bh$births_wt <- bh$wmweight
bh$HH6[bh$HH6=="Nomadic"] <- "Rural"


gfr <- merge(stats::aggregate(births_wt ~ HH7, bh, sum),
             stats::aggregate(pys_wt ~ HH7, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)
gfr
sum(gfr$births_wt) / sum(gfr$pys_wt)

subset(subp_gfr, country == "Somalia")
subp_gfr[subp_gfr$country == "Somalia" & subp_gfr$eppregion == "Puntland", c("survyear", "survtype", "gfr")] <- c(2006, "MICS", 0.204)
subp_gfr[subp_gfr$country == "Somalia" & subp_gfr$eppregion == "Somaliland", c("survyear", "survtype", "gfr")] <- c(2006, "MICS", 0.168)
subp_gfr[subp_gfr$country == "Somalia" & subp_gfr$eppregion == "South-Central", c("survyear", "survtype", "gfr")] <- c(2006, "MICS", 0.236)


## Central African Republic: MICS 2006

wm <- foreign::read.spss("~/Documents/Data/MICS/Central African Republic_MICS4_Datasets/Central African Republic MICS 2010 SPSS Datasets/wm.sav")
wm <- data.frame(wm)

wm$start <- pmax(wm$WDOI-36, wm$WDOB+15*12)
wm$end <- pmin(wm$WDOI, wm$WDOB+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$WMWEIGHT
wm <- subset(wm, pys_wt > 0)

bh <- foreign::read.spss("~/Documents/Data/MICS/Central African Republic_MICS4_Datasets/Central African Republic MICS 2010 SPSS Datasets/bh.sav")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("HH1", "HH2", "LN", "start", "end")], by=c("HH1", "HH2", "LN"))
bh <- subset(bh, IMPBH4C > start & IMPBH4C <= end)
bh$births_wt <- bh$wmweight
bh$HH6[bh$HH6=="Nomadic"] <- "Rural"


gfr <- merge(stats::aggregate(births_wt ~ HH7, bh, sum),
             stats::aggregate(pys_wt ~ HH7, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)
gfr
sum(gfr$births_wt) / sum(gfr$pys_wt)

subset(subp_gfr, country == "Central African Republic")
subp_gfr[subp_gfr$country == "Central African Republic" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.117)
subp_gfr[subp_gfr$country == "Central African Republic" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.106)
subp_gfr[subp_gfr$country == "Central African Republic" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2014, "MICS", 0.122)




## Botswana: Botswana Family Health Survey 2007-08

wm <- foreign::read.spss("~/Documents/Data/Botswana/Famliy Health Survey IV 2007-08/bfhs-2007-2008-wm-v1.sav")
wm <- data.frame(wm)

wm$start <- pmax(wm$cmcdoiw-36, wm$wdob+15*12)
wm$end <- pmin(wm$cmcdoiw, wm$wdob+45*12)
wm$pys <- (wm$end - wm$start)/12
wm$pys_wt <- wm$pys*wm$wmweight
wm <- subset(wm, pys_wt > 0)

bh <- foreign::read.spss("~/Documents/Data/Botswana/Famliy Health Survey IV 2007-08/bfhs-2007-2008-bh-v1.sav")
bh.variable.labels <- attr(bh, "variable.labels")
bh.label.tabel <- attr(bh, "label.table")
bh <- data.frame(bh)

bh <- merge(bh, wm[c("EA_ID", "HH_NO", "FE_LINE_NO", "start", "end")], by=c("EA_ID", "HH_NO", "FE_LINE_NO"))
bh <- subset(bh, ccdob > start & ccdob <= end)
bh$births_wt <- bh$wmweight

bh$residence[bh$residence == "City/Town"] <- "Urban Village"
wm$residence[wm$residence == "City/Town"] <- "Urban Village"

gfr <- merge(stats::aggregate(births_wt ~ residence, bh, sum),
             stats::aggregate(pys_wt ~ residence, wm, sum))
gfr$gfr <- with(gfr, births_wt / pys_wt)
gfr
sum(gfr$births_wt) / sum(gfr$pys_wt)


subset(subp_gfr, country == "Botswana")
subp_gfr[subp_gfr$country == "Botswana" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2008, "MICS", 0.103)
subp_gfr[subp_gfr$country == "Botswana" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2008, "MICS", 0.088)
subp_gfr[subp_gfr$country == "Botswana" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2008, "MICS", 0.130)



## Mauritius: WPP 2017 (//esa.un.org/unpd/wpp/)
## Use same value for Urban and Rural

subset(subp_gfr, country == "Mauritius")
subp_gfr[subp_gfr$country == "Mauritius" & subp_gfr$eppregion == "N", c("survyear", "survtype", "gfr")] <- c(2015, "WPP", 0.051)
subp_gfr[subp_gfr$country == "Mauritius" & subp_gfr$eppregion == "U", c("survyear", "survtype", "gfr")] <- c(2015, "WPP", 0.051)
subp_gfr[subp_gfr$country == "Mauritius" & subp_gfr$eppregion == "R", c("survyear", "survtype", "gfr")] <- c(2015, "WPP", 0.051)



#################################
####  Save subp_gfr dataset  ####
#################################

subp_gfr$gfr <- as.numeric(subp_gfr$gfr)
subp_gfr$survyear <- as.integer(subp_gfr$survyear)

devtools::use_data(subp, subp_gfr, internal = TRUE, overwrite=TRUE)
write.csv(subp_gfr, file="~/Downloads/subp_gfr.csv", row.names=FALSE)
