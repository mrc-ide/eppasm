get_dp_version <- function(dp){
  dp.vers <- dp[2,1] # <General 3>: 2013, 2014 Spectrum files; <General5>: 2015 Spectrum files
  if(!dp.vers %in% c("<General 3>", "<General5>"))
    if(dp.vers == "<FirstYear MV>")
      dp.vers <- "Spectrum2016"
    else if(dp.vers == "<FirstYear MV2>")
      dp.vers <- "Spectrum2017"
    else
      stop("Spectrum DP file version not recognized. Package probably needs to be updated to most recent Spectrum version.")
  return(dp.vers)
}

read_dp <- function(pjnz){
  dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)
  return(dp)
}

read_pjn <- function(pjnz){
  dpfile <- grep(".PJN$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)
  return(dp)
}


read_region <- function(pjnz){
  pjn <- read_pjn(pjnz)
  region <- pjn[which(pjn[,1] == "<Projection Parameters - Subnational Region Name2>")+2, 4]
  if(region == "")
    return(NULL)
  else
    return(region)
}

read_country <- function(pjnz){
  pjn <- read_pjn(pjnz)
  cc <- as.integer(pjn[which(pjn[,1] == "<Projection Parameters>")+2, 4])
  return(with(spectrum5_countrylist, Country[Code == cc]))
}  

###################################################
####  function to read HIV projection outputs  ####
###################################################

read_hivproj_output <- function(pjnz, single.age=TRUE){

  ## read .DP file
  dp <- read_dp(pjnz)

  dp.vers <- get_dp_version(dp)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}

  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  if(dp.vers %in% c("<General 3>", "<General5>")){
    version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
    validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
    validversion <- dp[which(dp[,1] == "<ValidVers>")+2,4]
  } else if(dp.vers == "Spectrum2016") {
    version <- as.numeric(dpsub("<VersionNum MV>", 3, 4))
    validdate <- dpsub("<ValidDate MV>",2,3)
    validversion <- dpsub("<ValidVers MV>",2,4)
  } else if(dp.vers == "Spectrum2017") {
    version <- as.numeric(dpsub("<VersionNum MV2>", 3, 4))
    validdate <- dpsub("<ValidDate MV>",2,3)
    validversion <- dpsub("<ValidVers MV>",2,4)
  }


  ## projection parameters
  if(dp.vers %in% c("<General 3>", "<General5>")){
    yr_start <- as.integer(dp[which(dp[,2] == "First year")+1,4])
    yr_end <- as.integer(dp[which(dp[,2] == "Final year")+1,4])
  } else if(dp.vers == "Spectrum2016"){
    yr_start <- as.integer(dpsub("<FirstYear MV>",3,4))
    yr_end <- as.integer(dpsub("<FinalYear MV>",3,4))
  } else if(dp.vers == "Spectrum2017"){
    yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
    yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
    t0 <- as.numeric(dpsub("<FirstYearOfEpidemic MV>",2,4))
  }
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1

  agegr.lab <- c(paste(0:15*5, 1:16*5, sep="-"), "80+")

  ## Number HIV+
  if(dp.vers %in% c("<General 3>", "<General5>")){
    hivnum.age <- sapply(dpsub("HIV", 4:54, timedat.idx, 2), as.numeric)
    rownames(hivnum.age) <- dpsub("HIV", 4:54, 3, 2)
  } else {
    hivnum.age <- sapply(dpsub("<HIV MV>", 5:55, timedat.idx), as.numeric)
    rownames(hivnum.age) <- dpsub("<HIV MV>", 5:55, 3)
  }
  hivnum.m <- hivnum.age[grep("Sex=1", rownames(hivnum.age)),]
  hivnum.f <- hivnum.age[grep("Sex=2", rownames(hivnum.age)),]
  dimnames(hivnum.m) <- dimnames(hivnum.f) <- list(agegr.lab, proj.years)

  ## Number new infections
  if(dp.vers %in% c("<General 3>", "<General5>")){
    newinf.age <- sapply(dpsub("New Infections", 4:54, timedat.idx, 2), as.numeric)
    rownames(newinf.age) <- dpsub("HIV", 4:54, 3, 2)
  } else {
    newinf.age <- sapply(dpsub("<NewInfections MV>", 5:55, timedat.idx), as.numeric)
    rownames(newinf.age) <- dpsub("<NewInfections MV>", 5:55, 3)
  }
  newinf.m <- newinf.age[grep("Sex=1", rownames(newinf.age)),]
  newinf.f <- newinf.age[grep("Sex=2", rownames(newinf.age)),]
  dimnames(newinf.m) <- dimnames(newinf.f) <- list(agegr.lab, proj.years)

  ## Number of infant (age 0) new infections

  ## Total population size
  if(dp.vers == "<General 3>"){
    totpop.tidx <- which(dp[,2] == "Total Population")
    totpop.m <- sapply(dp[totpop.tidx + 1:17*7 + 6, timedat.idx], as.numeric)
    totpop.f <- sapply(dp[totpop.tidx + 1:17*7 + 8, timedat.idx], as.numeric)
  } else if(dp.vers == "<General5>"){
    ## totpop.tidx <- which(dp[,1] == "<BigPop3>")
    totpop.m.tidx <- which(dp[,3] == "Males, Total, Age 0")
    totpop.m <- sapply(lapply(dp[totpop.m.tidx + 1:81 - 1, timedat.idx], as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f.tidx <- which(dp[,3] == "Females, Total, Age 0")
    totpop.f <- sapply(lapply(dp[totpop.f.tidx + 1:81 - 1, timedat.idx], as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
  } else if(exists_dptag("<BigPop MV>")){
    totpop.m <- sapply(lapply(dpsub("<BigPop MV>", 3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f <- sapply(lapply(dpsub("<BigPop MV>", 81+3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
  } else if(exists_dptag("<BigPop MV2>")){
    totpop.m <- sapply(lapply(dpsub("<BigPop MV2>", 3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f <- sapply(lapply(dpsub("<BigPop MV2>", 243+3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
  } else if(exists_dptag("<BigPop MV3>")){
    totpop.m <- sapply(lapply(dpsub("<BigPop MV3>", 3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f <- sapply(lapply(dpsub("<BigPop MV3>", 84:164, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
  }
  dimnames(totpop.m) <- dimnames(totpop.f) <- list(agegr.lab, proj.years)

  ## ART need
  if(dp.vers %in% c("<General 3>", "<General5>")){
    artneed.m <- sapply(dpsub("Need FL", 5+0:16*3, timedat.idx, 2), as.numeric)
    artneed.f <- sapply(dpsub("Need FL", 6+0:16*3, timedat.idx, 2), as.numeric)
  } else {
    artneed.m <- sapply(dpsub("<NeedART MV>", 6+0:16*3, timedat.idx), as.numeric)
    artneed.f <- sapply(dpsub("<NeedART MV>", 7+0:16*3, timedat.idx), as.numeric)
  }
  dimnames(artneed.m) <- dimnames(artneed.f) <- list(agegr.lab, proj.years)

  ## On ART
  if(dp.vers %in% c("<General 3>", "<General5>")){
    artnum.m <- sapply(dpsub("On FL", 5+0:16*3, timedat.idx, 2), as.numeric)
    artnum.f <- sapply(dpsub("On FL", 6+0:16*3, timedat.idx, 2), as.numeric)
  } else {
    artnum.m <- sapply(dpsub("<OnART MV>", 6+0:16*3, timedat.idx), as.numeric)
    artnum.f <- sapply(dpsub("<OnART MV>", 7+0:16*3, timedat.idx), as.numeric)
  }
  dimnames(artnum.m) <- dimnames(artnum.f) <- list(agegr.lab, proj.years)
  
  ## number non-aids deaths
  if(dp.vers %in% c("<General 3>", "<General5>")){
    natdeaths.m <- sapply(dpsub("Deaths - Male", 2:18, timedat.idx, 2), as.numeric)
    natdeaths.f <- sapply(dpsub("Deaths - Female", 2:18, timedat.idx, 2), as.numeric)
  } else if(dp.vers == "Spectrum2016"){
    natdeaths.m <- sapply(dpsub("<Deaths MV>", 5:21, timedat.idx), as.numeric)
    natdeaths.f <- sapply(dpsub("<Deaths MV>", 24:40, timedat.idx), as.numeric)
  } else if(dp.vers == "Spectrum2017"){
    natdeaths.m <- sapply(dpsub("<Deaths MV2>", 4:20, timedat.idx), as.numeric)
    natdeaths.f <- sapply(dpsub("<Deaths MV2>", 22:38, timedat.idx), as.numeric)
  }
  dimnames(natdeaths.m) <- dimnames(natdeaths.f) <- list(agegr.lab, proj.years)

  ## number AIDS deaths
  if(dp.vers %in% c("<General 3>", "<General5>")){
    aidsdeaths.m <- sapply(dpsub("AIDS deaths - Male", 2:18, timedat.idx, 2), as.numeric)
    aidsdeaths.f <- sapply(dpsub("AIDS deaths - Female", 2:18, timedat.idx, 2), as.numeric)
  } else if(dp.vers == "Spectrum2016"){
    aidsdeaths.m <- sapply(dpsub("<AIDSDeaths MV>", 3:19, timedat.idx), as.numeric)
    aidsdeaths.f <- sapply(dpsub("<AIDSDeaths MV>", 22:38, timedat.idx), as.numeric)
  } else if(dp.vers == "Spectrum2017"){
    aidsdeaths.m <- sapply(dpsub("<AIDSDeaths MV2>", 4:20, timedat.idx), as.numeric)
    aidsdeaths.f <- sapply(dpsub("<AIDSDeaths MV2>", 22:38, timedat.idx), as.numeric)
  }
  dimnames(aidsdeaths.m) <- dimnames(aidsdeaths.f) <- list(agegr.lab, proj.years)

  specres <- list("totpop.m" = totpop.m,
                  "totpop.f" = totpop.f,
                  "hivnum.m" = hivnum.m,
                  "hivnum.f" = hivnum.f,
                  "artnum.m" = artnum.m,
                  "artnum.f" = artnum.f,
                  "newinf.m" = newinf.m,
                  "newinf.f" = newinf.f,
                  "natdeaths.m" = natdeaths.m,
                  "natdeaths.f" = natdeaths.f,
                  "aidsdeaths.m" = aidsdeaths.m,
                  "aidsdeaths.f" = aidsdeaths.f)

  if(single.age){
    if(exists_dptag("<BigPop3>"))
      totpop <- sapply(dpsub("<BigPop3>", 2:163, timedat.idx), as.numeric)
    else if(exists_dptag("<BigPop MV>"))
      totpop <- sapply(dpsub("<BigPop MV>", 3:164, timedat.idx), as.numeric)
    else if(exists_dptag("<BigPop MV2>"))
      totpop <- sapply(dpsub("<BigPop MV2>", c(3+0:80, 246+0:80), timedat.idx), as.numeric)
    else if(exists_dptag("<BigPop MV3>"))
      totpop <- sapply(dpsub("<BigPop MV3>", 3:164, timedat.idx), as.numeric)
    totpop <- array(totpop, c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
                      
    if(exists_dptag("<HIVBySingleAge MV>"))
      hivpop <- array(sapply(dpsub("<HIVBySingleAge MV>", c(3:83, 85:165), timedat.idx), as.numeric),
                      c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
    else if(exists_dptag("<HIVBySingleAge MV2>"))
      hivpop <- sapply(dpsub("<HIVBySingleAge MV2>", 3:164, timedat.idx), as.numeric)
    else
      hivpop <- NA
    hivpop <- array(hivpop, c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
    
    if(exists_dptag("<DeathsByAge MV>"))
      natdeaths <- sapply(dpsub("<DeathsByAge MV>", c(4:84, 86:166), timedat.idx), as.numeric)
    else if(exists_dptag("<DeathsByAge MV2>"))
      natdeaths <- sapply(dpsub("<DeathsByAge MV2>", 3:164, timedat.idx), as.numeric)
    else
      natdeaths <- NA
    natdeaths <- array(natdeaths, c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
    
    if(exists_dptag("<AidsDeathsByAge MV>"))
      hivdeaths <- sapply(dpsub("<AidsDeathsByAge MV>", c(4:84, 86:166), timedat.idx), as.numeric)
    else if(exists_dptag("<AidsDeathsByAge MV2>"))
      hivdeaths <- sapply(dpsub("<AidsDeathsByAge MV2>", 3:164, timedat.idx), as.numeric)
    else
      hivdeaths <- NA
    hivdeaths <- array(hivdeaths, c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
    
    specres[c("totpop", "hivpop", "natdeaths", "hivdeaths")] <- list(totpop, hivpop, natdeaths, hivdeaths)
  }

  specres$births <- setNames(as.numeric(dpsub("<Births MV>", 2, timedat.idx)), proj.years)
  specres$hivpregwomen <- setNames(as.numeric(dpsub("<ChildNeedPMTCT MV>", 2, timedat.idx)), proj.years)
  specres$hivpregwomen <- setNames(as.numeric(dpsub("<ChildNeedPMTCT MV>", 2, timedat.idx)), proj.years)
  specres$receivepmtct <- setNames(as.numeric(dpsub("<ChildOnPMTCT MV>", 2, timedat.idx)), proj.years)
  
  class(specres) <- "specres"
  attr(specres, "country") <- read_country(pjnz)
  attr(specres, "region") <- read_region(pjnz)

  return(specres)
}


######################################################
####  function to read HIV projection parameters  ####
######################################################

read_hivproj_param <- function(pjnz, use_ep5=FALSE){

  ## read .DP file

  if(use_ep5)
    dpfile <- grep(".ep5$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  else
    dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  if(use_ep5)
    dp.vers <- "Spectrum2017"
  else
    dp.vers <- get_dp_version(dp)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}

  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  if(dp.vers %in% c("<General 3>", "<General5>")){
    version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
    validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
    validversion <- as.numeric(dp[which(dp[,1] == "<ValidVers>")+2,4])
  } else if(dp.vers == "Spectrum2016") {
    version <- as.numeric(dpsub("<VersionNum MV>", 2, 4))
    validdate <- dpsub("<ValidDate MV>",2,3)
    validversion <- dpsub("<ValidVers MV>",2,4)
  } else if(dp.vers == "Spectrum2017") {
    version <- as.numeric(dpsub("<VersionNum MV2>", 3, 4))
    validdate <- dpsub("<ValidDate MV>",2,3)
    validversion <- dpsub("<ValidVers MV>",2,4)
  }


  ## find tag indexes (tidx)
  if(dp.vers == "<General 3>"){
    aids5.tidx <- which(dp[,1] == "<AIDS5>")
  } else if(dp.vers == "<General5>"){
    hivtfr.tidx <- which(dp[,1] == "<HIVTFR2>")
    hivsexrat.tidx <- which(dp[,1] == "<HIVSexRatio>")
    hivagedist.tidx <- which(dp[,1] == "<HIVDistribution2>")
  }

  if(dp.vers %in% c("<General 3>", "<General5>")){
    nathist.tidx <- which(dp[,1] == "<AdultTransParam2>")
    cd4initdist.tidx <- which(dp[,1] == "<DistNewInfectionsCD4>")
    infectreduc.tidx <- which(dp[,1] == "<InfectReduc>")
    adult.artnumperc.tidx <- which(dp[,1] == "<HAARTBySexPerNum>")
    adult.art.tidx <- which(dp[,1] == "<HAARTBySex>")
    adult.arteligthresh.tidx <- which(dp[,1] == "<CD4ThreshHoldAdults>")
    specpopelig.tidx <- which(dp[,1] == "<PopsEligTreat1>")
  }

  ## state space dimensions
  NG <- 2
  AG <- 17
  DS <- 7
  TS <- 3
  
  ## projection parameters
  if(dp.vers %in% c("<General 3>", "<General5>")){
    yr_start <- as.integer(dp[which(dp[,2] == "First year")+1,4])
    yr_end <- as.integer(dp[which(dp[,2] == "Final year")+1,4])
  } else if(dp.vers == "Spectrum2016"){
    yr_start <- as.integer(dpsub("<FirstYear MV>",3,4))
    yr_end <- as.integer(dpsub("<FinalYear MV>",3,4))
  } else if(dp.vers == "Spectrum2017"){
    yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
    yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  }
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
  
  ## scalar paramters
  if(dp.vers %in% c("<General 3>", "<General5>")){
     relinfectART <- 1.0 - as.numeric(dp[infectreduc.tidx+1, 4])
  } else if(dp.vers %in% c("Spectrum2016", "Spectrum2017")){
    relinfectART <- 1.0 - as.numeric(dpsub("<AdultInfectReduc MV>",2,4))
  }

  if(dp.vers == "<General 3>"){
    fert_rat <- as.numeric(dp[which(dp[,1] == "<AIDS5>")+185, 4+0:6])
    fert_rat <- array(rep(fert_rat, length(proj.years)), c(7, length(proj.years)))
    dimnames(fert_rat) <- list(agegr=seq(15, 45, 5), year=proj.years)
  } else if(dp.vers == "<General5>") {
    fert_rat <- sapply(dp[hivtfr.tidx+2:8, 3+seq_along(proj.years)], as.numeric)
    dimnames(fert_rat) <- list(agegr=seq(15, 45, 5), year=proj.years)
  } else if(dp.vers == "Spectrum2016") {
    fert_rat <- sapply(dpsub("<HIVTFR MV>", 2:8, timedat.idx), as.numeric)
    dimnames(fert_rat) <- list(agegr=seq(15, 45, 5), year=proj.years)
  } else if(exists_dptag("<HIVTFR MV2>")) {
    fert_rat <- sapply(dpsub("<HIVTFR MV2>", 2:7, timedat.idx), as.numeric)
    dimnames(fert_rat) <- list(agegr=c(15, 18, seq(20, 35, 5)), year=proj.years)  # this version of Spectrum stratified fertility reduction by 15-17, 18-19, 20-24, ...
  } else if(exists_dptag("<HIVTFR MV3>")){
    fert_rat <- sapply(dpsub("<HIVTFR MV3>", 2:8, timedat.idx), as.numeric)
    dimnames(fert_rat) <- list(agegr=seq(15, 45, 5), year=proj.years)
  }

  if(dp.vers == "Spectrum2017")
    cd4fert_rat <- as.numeric(dpsub("<FertCD4Discount MV>", 2, 4+1:DS))
  else
    cd4fert_rat <- rep(1.0, DS)

  if(exists_dptag("<RatioWomenOnART MV>"))
    frr_art6mos <- as.numeric(dpsub("<RatioWomenOnART MV>", 2, 4))
  else
    frr_art6mos <- 1.0

  ## sex/age-specific incidence ratios (time varying)
  incrr_age <- array(NA, c(AG, NG, length(proj.years)), list(0:(AG-1)*5, c("Male", "Female"), proj.years))
  if(dp.vers == "<General 3>"){
    incrr_sex <- setNames(as.numeric(dp[aids5.tidx+181,timedat.idx]), proj.years) # !!! Not sure what aids5.tidx+183 (Ratio of female to male prevalence)
    incrr_age[,"Male",] <- sapply(dp[aids5.tidx+281:297,timedat.idx], as.numeric)
    incrr_age[,"Female",] <- sapply(dp[aids5.tidx+299:315,timedat.idx], as.numeric)
  } else if(dp.vers == "<General5>"){
    incrr_sex <- setNames(as.numeric(dp[hivsexrat.tidx+2, timedat.idx]), proj.years)
    incrr_age[,"Male",] <- sapply(dp[hivagedist.tidx+3:19,timedat.idx], as.numeric)
    incrr_age[,"Female",] <- sapply(dp[hivagedist.tidx+21:37,timedat.idx], as.numeric)
  } else if(dp.vers == "Spectrum2016") {
    incrr_sex <- setNames(as.numeric(dpsub("<HIVSexRatio MV>", 2, timedat.idx)), proj.years)
    incrr_age[,"Male",] <- sapply(dpsub("<DistOfHIV MV>", 4:20, timedat.idx), as.numeric)
    incrr_age[,"Female",] <- sapply(dpsub("<DistOfHIV MV>", 22:38, timedat.idx), as.numeric)
  } else if(dp.vers == "Spectrum2017") {
    incrr_sex <- setNames(as.numeric(dpsub("<HIVSexRatio MV>", 2, timedat.idx)), proj.years)
    incrr_age[,"Male",] <- sapply(dpsub("<DistOfHIV MV2>", 3:19, timedat.idx), as.numeric)
    incrr_age[,"Female",] <- sapply(dpsub("<DistOfHIV MV2>", 20:36, timedat.idx), as.numeric)
  }


  ## hiv natural history
  cd4_initdist <- array(NA, c(DS, 4, NG), list(cd4stage=1:DS, agecat=c("15-24", "25-34", "35-44", "45+"), sex=c("Male", "Female")))
  cd4_prog <- array(NA, c(DS-1, 4, NG), list(cd4stage=1:(DS-1), agecat=c("15-24", "25-34", "35-44", "45+"), sex=c("Male", "Female")))
  cd4_mort <- array(NA, c(DS, 4, NG), list(cd4stage=1:DS, agecat=c("15-24", "25-34", "35-44", "45+"), sex=c("Male", "Female")))
  art_mort <- array(NA, c(TS, DS, 4, NG), list(artdur=c("ART0MOS", "ART6MOS", "ART1YR"), cd4stage=1:DS, agecat=c("15-24", "25-34", "35-44", "45+"), sex=c("Male", "Female")))

  if(dp.vers %in% c("<General 3>", "<General5>")){
    cd4_initdist[,,"Male"] <- array(as.numeric(dp[cd4initdist.tidx+2, 4:31])/100, c(DS, 4))
    cd4_initdist[,,"Female"] <- array(as.numeric(dp[cd4initdist.tidx+3, 4:31])/100, c(DS, 4))

    cd4_prog[,,"Male"] <- array(1/as.numeric(dp[nathist.tidx+4, 4:27]), c(DS-1, 4))
    cd4_prog[,,"Female"] <- array(1/as.numeric(dp[nathist.tidx+5, 4:27]), c(DS-1, 4))

    cd4_mort[,,"Male"] <- array(as.numeric(dp[nathist.tidx+7, 4:31]), c(DS, 4))
    cd4_mort[,,"Female"] <- array(as.numeric(dp[nathist.tidx+8, 4:31]), c(DS, 4))

    art_mort[1,,,"Male"] <- array(as.numeric(dp[nathist.tidx+10, 4:31]), c(DS, 4))
    art_mort[1,,,"Female"] <- array(as.numeric(dp[nathist.tidx+11, 4:31]), c(DS, 4))
    art_mort[2,,,"Male"] <- array(as.numeric(dp[nathist.tidx+13, 4:31]), c(DS, 4))
    art_mort[2,,,"Female"] <- array(as.numeric(dp[nathist.tidx+14, 4:31]), c(DS, 4))
    art_mort[3,,,"Male"] <- array(as.numeric(dp[nathist.tidx+16, 4:31]), c(DS, 4))
    art_mort[3,,,"Female"] <- array(as.numeric(dp[nathist.tidx+17, 4:31]), c(DS, 4))
    
  } else if(dp.vers %in% c("Spectrum2016", "Spectrum2017")) {
    cd4_initdist[,,"Male"] <- array(as.numeric(dpsub("<AdultDistNewInfectionsCD4 MV>", 3, 4:31))/100, c(DS, 4))
    cd4_initdist[,,"Female"] <- array(as.numeric(dpsub("<AdultDistNewInfectionsCD4 MV>", 4, 4:31))/100, c(DS, 4))

    ## Note: CD4 progression array has DS values, but should only be DS-1. Not sure what the last one is.
    cd4_prog[,,"Male"] <- array(as.numeric(dpsub("<AdultAnnRateProgressLowerCD4 MV>", 3, 4:31)), c(DS, 4))[1:(DS-1),]
    cd4_prog[,,"Female"] <- array(as.numeric(dpsub("<AdultAnnRateProgressLowerCD4 MV>", 4, 4:31)), c(DS, 4))[1:(DS-1),]

    cd4_mort[,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4NoART MV>", 3, 4:31)), c(DS, 4))
    cd4_mort[,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4NoART MV>", 4, 4:31)), c(DS, 4))
    
    if(dp.vers == "Spectrum2016"){
      art_mort[1,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV>", 3, 4:31)), c(DS, 4))
      art_mort[1,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV>", 4, 4:31)), c(DS, 4))
      art_mort[2,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV>", 3, 4:31)), c(DS, 4))
      art_mort[2,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV>", 4, 4:31)), c(DS, 4))
      art_mort[3,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV>", 3, 4:31)), c(DS, 4))
      art_mort[3,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV>", 4, 4:31)), c(DS, 4))
    } else if(dp.vers == "Spectrum2017") {
      art_mort[1,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[1,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV2>", 3, 4:31)), c(DS, 4))
      art_mort[2,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[2,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV2>", 3, 4:31)), c(DS, 4))
      art_mort[3,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[3,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV2>", 3, 4:31)), c(DS, 4))
    }
  }

  
  ## program parameters
  if(dp.vers %in% c("<General 3>", "<General5>")){
    art15plus_numperc <- sapply(dp[adult.artnumperc.tidx+3:4, timedat.idx], as.numeric)
    art15plus_num <- sapply(dp[adult.art.tidx+3:4, timedat.idx], as.numeric)
    art15plus_eligthresh <- setNames(as.numeric(dp[adult.arteligthresh.tidx+2, timedat.idx]), proj.years)
    artelig_specpop <- setNames(dp[specpopelig.tidx+1:7,2:6], c("description", "pop", "elig", "percent", "year"))
  } else if(dp.vers %in% c("Spectrum2016", "Spectrum2017")) {
    art15plus_numperc <- sapply(dpsub("<HAARTBySexPerNum MV>", 4:5, timedat.idx), as.numeric)
    art15plus_num <- sapply(dpsub("<HAARTBySex MV>", 4:5, timedat.idx), as.numeric)
    art15plus_eligthresh <- setNames(as.numeric(dpsub("<CD4ThreshHoldAdults MV>", 2, timedat.idx)), proj.years)
    artelig_specpop <- setNames(dpsub("<PopsEligTreat MV>", 3:9, 2:6), c("description", "pop", "elig", "percent", "year"))
  }
    
  dimnames(art15plus_numperc) <- list(sex=c("Male", "Female"), year=proj.years)
  dimnames(art15plus_num) <- list(sex=c("Male", "Female"), year=proj.years)

  artelig_specpop$pop <- c("PW", "TBHIV", "DC", "FSW", "MSM", "IDU", "OTHER")
  artelig_specpop$elig <- as.logical(as.integer(artelig_specpop$elig))
  artelig_specpop$percent <- as.numeric(artelig_specpop$percent)/100
  artelig_specpop$year <- as.integer(artelig_specpop$year)
  artelig_specpop$idx <- match(as.integer(artelig_specpop$year), proj.years)
  rownames(artelig_specpop) <- artelig_specpop$pop

  if(dp.vers %in% c("Spectrum2016", "Spectrum2017"))
    median_cd4init <- sapply(dpsub("<MedCD4CountInit MV>", 2, timedat.idx), as.numeric)
  else
    median_cd4init <- rep(0, length(timedat.idx))
  names(median_cd4init) <- proj.years
  
  if(dp.vers %in% c("Spectrum2016", "Spectrum2017"))
    art_dropout <- sapply(dpsub("<PercLostFollowup MV>", 2, timedat.idx), as.numeric)
  else
    art_dropout <- rep(0, length(timedat.idx))
  names(art_dropout) <- proj.years


  ## vertical transmission and paediatric survival

  verttrans <- setNames(sapply(dpsub("<PerinatalTransmission MV>", 2, timedat.idx), as.numeric)/100, proj.years)

  if(exists_dptag("<HIVBySingleAge MV>"))
    hivpop <- array(sapply(dpsub("<HIVBySingleAge MV>", c(3:83, 85:165), timedat.idx), as.numeric),
                    c(81, NG, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else if(exists_dptag("<HIVBySingleAge MV2>"))
    hivpop <- array(sapply(dpsub("<HIVBySingleAge MV2>", 3:164, timedat.idx), as.numeric),
                    c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else
    hivpop <- NULL

  if(exists_dptag("<AidsDeathsByAge MV>"))
    hivdeaths <- array(sapply(dpsub("<AidsDeathsByAge MV>", c(4:84, 86:166), timedat.idx), as.numeric),
                       c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else if(exists_dptag("<AidsDeathsByAge MV2>"))
    hivdeaths <- array(sapply(dpsub("<AidsDeathsByAge MV2>", 3:164, timedat.idx), as.numeric),
                       c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else
    hivpop <- NULL

  ## distribution of age 14 population
  PAED_DS <- 6 # number of paediatric stages of infection
  if(exists_dptag("<ChAged14ByCD4Cat MV>")){
    age14hivpop <- sapply(dpsub("<ChAged14ByCD4Cat MV>", 1+1:(NG*PAED_DS*(4+TS)), timedat.idx), as.numeric)
    age14hivpop <- array(age14hivpop, c(4+TS, PAED_DS, NG, length(proj.years)),
                         list(ARTstage=c("PERINAT", "BF0MOS", "BF6MOS", "BF1YR", "ART0MOS", "ART6MOS", "ART1YR"),
                              CD4cat=c("CD4_1000", "CD4_750", "CD4_500", "CD4_350", "CD4_200", "CD4_0"),
                              Sex=c("Male", "Female"), Year=proj.years))
  } else {
    
    ## Approximate for versions of Spectrum < 5.63
    specres <- read_hivproj_output(pjnz)
    hivpop14 <- specres$hivpop["14",,]

    ## Assume ART coverage for age 10-14 age group
    artcov14 <- rbind(Male = specres$artnum.m["10-15",]/specres$hivnum.m["10-15",],
                      Female = specres$artnum.f["10-15",]/specres$hivnum.f["10-15",])
    artcov14[is.na(artcov14)] <- 0
                      
    noart_cd4dist <- c(0.01, 0.04, 0.12, 0.22, 0.26, 0.35) # approximation for pre-ART period

    age14hivpop <- array(0, c(4+TS, PAED_DS, NG, length(proj.years)),
                         list(ARTstage=c("PERINAT", "BF0MOS", "BF6MOS", "BF1YR", "ART0MOS", "ART6MOS", "ART1YR"),
                              CD4cat=c("CD4_1000", "CD4_750", "CD4_500", "CD4_350", "CD4_200", "CD4_0"),
                              Sex=c("Male", "Female"), Year=proj.years))

    age14hivpop["PERINAT",,,] <- noart_cd4dist %o% (hivpop14 * (1 - artcov14))
    age14hivpop["ART1YR", "CD4_0",,] <- hivpop14 * artcov14
  }

  if(exists_dptag("<BigPop3>"))
    totpop <- sapply(dpsub("<BigPop3>", 2:163, timedat.idx), as.numeric)
  else if(exists_dptag("<BigPop MV>"))
    totpop <- sapply(dpsub("<BigPop MV>", 3:164, timedat.idx), as.numeric)
  else if(exists_dptag("<BigPop MV2>"))
    totpop <- sapply(dpsub("<BigPop MV2>", c(3+0:80, 246+0:80), timedat.idx), as.numeric)
  else if(exists_dptag("<BigPop MV3>"))
    totpop <- sapply(dpsub("<BigPop MV3>", 3:164, timedat.idx), as.numeric)
  else
    totpop <- NULL

  if(!is.null(totpop)){
    totpop <- array(totpop, c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
    age14totpop <- totpop["14",,]
  } else
    age14totpop <- NULL
  
  projp <- list("yr_start"=yr_start, "yr_end"=yr_end,
                "relinfectART"=relinfectART,
                "fert_rat"=fert_rat,
                "cd4fert_rat"=cd4fert_rat,
                "frr_art6mos"=frr_art6mos,
                "incrr_sex"=incrr_sex, "incrr_age"=incrr_age,
                "cd4_initdist"=cd4_initdist, "cd4_prog"=cd4_prog, "cd4_mort"=cd4_mort, "art_mort"=art_mort,
                "art15plus_numperc"=art15plus_numperc, "art15plus_num"=art15plus_num,
                "art15plus_eligthresh"=art15plus_eligthresh, "artelig_specpop"=artelig_specpop,
                "median_cd4init"=median_cd4init, "art_dropout"=art_dropout,
                "verttrans"=verttrans, "hivpop"=hivpop, "hivdeaths"=hivdeaths,
                "age14hivpop"=age14hivpop,
                "age14totpop"=age14totpop)
  class(projp) <- "projp"
  attr(projp, "version") <- version
  attr(projp, "validdate") <- validdate
  attr(projp, "validversion") <- validversion
  
  return(projp)
}


###################################################################
####  function to read UN Population Division projection file  ####
###################################################################


read_demog_param <- function(upd.file, age.intervals = 1){

  ## check age intervals and prepare age groups vector
  if(length(age.intervals) == 1){
    if(80 %% age.intervals != 0)
      stop("Invalid age interval interval")
    age.groups <- c(rep(seq(0, 79, age.intervals), each=age.intervals), 80)
    age.intervals <- c(rep(age.intervals, 80 / age.intervals), 1)
    ## } else if(length(age.groups) < 81 && age.groups[1] == 0 && tail(age.groups, 1) == 80){
    ##   stop("not defined yet") ## define thie one
  } else if(length(age.intervals) < 81){
    if(sum(age.intervals != 81))
      stop("Invalid vector of age intervals")
    age.groups <- rep(1:length(age.intervals), times=age.intervals)
  }

  ## Read and parse udp file
  upd <- read.csv(upd.file, header=FALSE, as.is=TRUE)

  wpp.version <- ifelse(any(upd[,1] ==  "RevisionYear=2015"), "wpp2015", "wpp2012")

  bp.tidx <- which(upd[,1] == "<basepop>")   # tag index [tidx]
  bp.eidx <- which(upd[,1] == "</basepop>")  # end tag index [eidx]

  lt.tidx <- which(upd[,1] == "<lfts>")
  lt.eidx <- which(upd[,1] == "</lfts>")

  pasfrs.tidx <- which(upd[,1] == "<pasfrs>")
  pasfrs.eidx <- which(upd[,1] == "</pasfrs>")

  migration.tidx <- which(upd[,1] == "<migration>")
  migration.eidx <- which(upd[,1] == "</migration>")

  tfr.tidx <- which(upd[,1] == "<tfr>")
  tfr.eidx <- which(upd[,1] == "</tfr>")

  srb.tidx <- which(upd[,1] == "<srb>")
  srb.eidx <- which(upd[,1] == "</srb>")

  bp <- setNames(data.frame(upd[(bp.tidx+2):(bp.eidx-1),1:4]), upd[bp.tidx+1,1:4])
  lt <- setNames(data.frame(upd[(lt.tidx+2):(lt.eidx-1),]), upd[lt.tidx+1,])
  pasfrs <- setNames(data.frame(upd[(pasfrs.tidx+2):(pasfrs.eidx-1),1:3]), upd[pasfrs.tidx+1,1:3])
  migration <- setNames(data.frame(upd[(migration.tidx+2):(migration.eidx-1),1:4]), upd[migration.tidx+1,1:4])
  tfr <- setNames(as.numeric(upd[(tfr.tidx+2):(tfr.eidx-1),2]), upd[(tfr.tidx+2):(tfr.eidx-1),1])
  srb <- setNames(as.numeric(upd[(srb.tidx+2):(srb.eidx-1),2]), upd[(srb.tidx+2):(srb.eidx-1),1])


  ## Aggregate into specified age groups

  ## population size
  basepop <- array(as.numeric(bp$value), c(length(unique(bp$age)), length(unique(bp$sex)), length(unique(bp$year))))
  dimnames(basepop) <- list(age=unique(bp$age), sex=c("Male", "Female"), year=unique(bp$year))
  basepop <- apply(basepop, 2:3, tapply, age.groups, sum)

  ## mx
  years <- unique(lt$year)
  nyears <- length(years)
  Sx <- as.numeric(lt$Sx[-(1:(2*nyears)*82-1)]) # 80+ age group given twice
  dim(Sx) <- c(81, 2, nyears)
  dimnames(Sx) <- list(age=0:80, sex=c("Male", "Female"), year=years)
  Sx <- apply(Sx, 2:3, tapply, age.groups, prod)
  mx <- -sweep(log(Sx), 1, age.intervals, "/")

  ## asfr
  asfd <- array(as.numeric(pasfrs$value), c(35, nyears))
  dimnames(asfd) <- list(age=15:49, year=years)
  asfr <- sweep(asfd, 2, tfr, "*")
  asfr <- apply(asfr, 2, tapply, age.groups[16:50], mean)

  asfd <- apply(asfd, 2, tapply, age.groups[16:50], sum)

  ## migration
  netmigr <- array(as.numeric(migration$value), c(81, 2, nyears))
  dimnames(netmigr) <- list(age=0:80, sex=c("Male", "Female"), year=years)
  netmigr <- apply(netmigr, 2:3, tapply, age.groups, sum)

  demp <- list("basepop"=basepop, "mx"=mx, "Sx"=Sx, "asfr"=asfr, "tfr"=tfr, "asfd"=asfd, "srb"=srb, "netmigr"=netmigr)
  class(demp) <- "demp"
  attr(demp, "version") <- wpp.version

  return(demp)
}


##########################################################
####  Read demographic inputs from Spectrum .DP file  ####
##########################################################

## Note: only parses Spectrum 2016 files, produces outputs by single-year age

read_specdp_demog_param <- function(pjnz, use_ep5=FALSE){

  if(use_ep5)
    dpfile <- grep(".ep5$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  else
    dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  if(use_ep5)
    dp.vers <- "Spectrum2017"
  else
    dp.vers <- get_dp_version(dp)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}
  
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  version <- paste("Spectrum", dp[which(dp[,1] == "<ValidVers MV>")+2, 4])

  ## projection parameters
  if(dp.vers == "Spectrum2016"){
    yr_start <- as.integer(dp[which(dp[,1] == "<FirstYear MV>")+3,4])
    yr_end <- as.integer(dp[which(dp[,1] == "<FinalYear MV>")+3,4])
  } else if(dp.vers == "Spectrum2017"){
    yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
    yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  } else
    stop(paste("Demographic inputs not available from Spectrum DP file (dp.vers =", dp.vers))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1


  ## population size
  if(exists_dptag("<BigPop MV>"))
    basepop <- array(sapply(dpsub("<BigPop MV>", 3:164, timedat.idx), as.numeric),
                    c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else if(exists_dptag("<BigPop MV2>"))
    basepop <- array(sapply(dpsub("<BigPop MV2>", c(3+0:80, 246+0:80), timedat.idx), as.numeric),
                     c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else if(exists_dptag("<BigPop MV3>"))
    basepop <- array(sapply(dpsub("<BigPop MV3>", 3+0:161, timedat.idx), as.numeric),
                     c(81, 2, length(proj.years)), list(0:80, c("Male", "Female"), proj.years))
  else
    stop("No recognized <BigPop MV[X]> tag, basepop not found.")
    
  
  ## mx
  if(dp.vers == "Spectrum2016"){
    sx.tidx <- which(dp[,1] == "<SurvRate MV>")
    Sx <- dp[sx.tidx+3+c(0:79,81, 83+0:79, 83+81), timedat.idx]
  } else if(dp.vers == "Spectrum2017")
    Sx <- dpsub("<SurvRate MV2>", 3+c(0:79, 81, 82+0:79, 82+81), timedat.idx)
  Sx <- array(as.numeric(unlist(Sx)), c(81, 2, length(proj.years)))
  dimnames(Sx) <- list(age=0:80, sex=c("Male", "Female"), year=proj.years)

  mx <- -log(Sx)

  ## asfr
  tfr.tidx <- which(dp[,1] == "<TFR MV>")
  asfd.tidx <- which(dp[,1] == "<ASFR MV>")

  tfr <- setNames(as.numeric(dp[tfr.tidx + 2, timedat.idx]), proj.years)
  asfd <- sapply(dp[asfd.tidx + 3:9, timedat.idx], as.numeric)/100
  asfd <- apply(asfd / 5, 2, rep, each=5)
  dimnames(asfd) <- list(age=15:49, year=proj.years)
  asfr <- sweep(asfd, 2, tfr, "*")

  births.tidx <- which(dp[,1] == "<Births MV>")
  births <- setNames(as.numeric(dp[births.tidx + 2, timedat.idx]), proj.years)

  ## srb
  srb.tidx <- which(dp[,1] == "<SexBirthRatio MV>")
  srb <- setNames(as.numeric(dp[srb.tidx + 2, timedat.idx]), proj.years)

  ## migration
  if(dp.vers == "Spectrum2016"){
    migrrate.tidx <- which(dp[,1] == "<MigrRate MV>")
    totnetmig <- sapply(dp[migrrate.tidx+c(5,8), timedat.idx], as.numeric)
  } else if(dp.vers == "Spectrum2017")
    totnetmig <- sapply(dpsub("<MigrRate MV2>", c(4, 6), timedat.idx), as.numeric)

  if(dp.vers == "Spectrum2016"){
    ## note: age=0 is empty in DP file, inputs start in age=1
    migaged.tidx <- which(dp[,1] == "<MigrAgeDist MV>")
    netmigagedist <- sapply(dp[migaged.tidx+6+c(1:17*2, 37+1:17*2), timedat.idx], as.numeric) / 100
  } else if(dp.vers == "Spectrum2017")
    netmigagedist <- sapply(dpsub("<MigrAgeDist MV2>", 2+1:34, timedat.idx), as.numeric) / 100
  netmigagedist <- array(c(netmigagedist), c(17, 2, length(proj.years)))
    
  netmigr <- sweep(netmigagedist, 2:3, totnetmig, "*")


  ## Beer's coefficients for disaggregating 5 year age groups into
  ## single-year age groups (from John Stover)
  Afirst <- rbind(c(0.3333, -0.1636, -0.0210,  0.0796, -0.0283),
                  c(0.2595, -0.0780,  0.0130,  0.0100, -0.0045),
                  c(0.1924,  0.0064,  0.0184, -0.0256,  0.0084),
                  c(0.1329,  0.0844,  0.0054, -0.0356,  0.0129),
                  c(0.0819,  0.1508, -0.0158, -0.0284,  0.0115))
  Asecond <- rbind(c( 0.0404,  0.2000, -0.0344, -0.0128,  0.0068),
                   c( 0.0093,  0.2268, -0.0402,  0.0028,  0.0013),
                   c(-0.0108,  0.2272, -0.0248,  0.0112, -0.0028),
                   c(-0.0198,  0.1992,  0.0172,  0.0072, -0.0038),
                   c(-0.0191,  0.1468,  0.0822, -0.0084, -0.0015))
  Amid <- rbind(c(-0.0117,  0.0804,  0.1570, -0.0284,  0.0027),
                c(-0.0020,  0.0160,  0.2200, -0.0400,  0.0060),
                c( 0.0050, -0.0280,  0.2460, -0.0280,  0.0050),
                c( 0.0060, -0.0400,  0.2200,  0.0160, -0.0020),
                c( 0.0027, -0.0284,  0.1570,  0.0804, -0.0117))
  Apenult <- rbind(c(-0.0015, -0.0084,  0.0822,  0.1468, -0.0191),
                   c(-0.0038,  0.0072,  0.0172,  0.1992, -0.0198),
                   c(-0.0028,  0.0112, -0.0248,  0.2272, -0.0108),
                   c( 0.0013,  0.0028, -0.0402,  0.2268,  0.0093),
                   c( 0.0068, -0.0128, -0.0344,  0.2000,  0.0404))
  Aultim <- rbind(c( 0.0115, -0.0284, -0.0158,  0.1508,  0.0819),
                  c( 0.0129, -0.0356,  0.0054,  0.0844,  0.1329),
                  c( 0.0084, -0.0256,  0.0184,  0.0064,  0.1924),
                  c(-0.0045,  0.0100,  0.0130, -0.0780,  0.2595),
                  c(-0.0283,  0.0796, -0.0210, -0.1636,  0.3333))

  A <- do.call(rbind,
               c(list(cbind(Afirst, matrix(0, 5, 12)),
                      cbind(Asecond, matrix(0, 5, 12))),
                 lapply(0:11, function(i) cbind(matrix(0, 5, i), Amid, matrix(0, 5, 12-i))),
                 list(cbind(matrix(0, 5, 11), Apenult, matrix(0, 5, 1)),
                      cbind(matrix(0, 5, 11), Aultim, matrix(0, 5, 1)),
                      c(rep(0, 16), 1))))

  netmigr <- apply(netmigr, 2:3, function(x) A %*% x)
  dimnames(netmigr) <- list(age=0:80, sex=c("Male", "Female"), year=proj.years)


  demp <- list("basepop"=basepop, "mx"=mx, "Sx"=Sx, "asfr"=asfr, "tfr"=tfr, "asfd"=asfd, "srb"=srb, "netmigr"=netmigr,
               "births"=births)
  class(demp) <- "demp"
  attr(demp, "version") <- version

  return(demp)
}



## Read percentage urban input from EPP XML file
#'
#' @param pjnz file path to Spectrum PJNZ file.
read_epp_perc_urban <- function(pjnz){

  xmlfile <- grep(".xml", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  con <- unz(pjnz, xmlfile)
  epp.xml <- scan(con, "character", sep="\n")
  close(con)
  
  if (!require("XML", quietly = TRUE))
    stop("read_epp_perc_urban() requires the package 'XML'. Please install it.", call. = FALSE)
  
  obj <- xmlTreeParse(epp.xml)
  r <- xmlRoot(obj)[[1]]

  yr_start <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetStartYear")]][[1]]))
  yr_end <- as.integer(xmlToList(r[[which(xmlSApply(r, xmlAttrs) == "worksetEndYear")]][[1]]))
  perc_urban.idx <- which(xmlSApply(r, xmlAttrs) == "currentUrbanPercent")
  if(length(perc_urban.idx) == 0){
    warning(paste0("EPP file does not contain Urban/Rural stratification:\n", pjnz))
    return(NULL)
  }
  perc_urban <- as.numeric(xmlSApply(r[[perc_urban.idx]][[1]], xmlSApply, xmlToList))

  return(setNames(perc_urban, yr_start:yr_end))
}
    
## Read epidemic start year from EPP XML file
#'
#' @param pjnz file path to Spectrum PJNZ file.
#' @return vector of epidemic start year for each EPP subregion with region names
read_epp_t0 <- function(pjnz){
  
  xmlfile <- grep(".xml", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  con <- unz(pjnz, xmlfile)
  epp.xml <- scan(con, "character", sep="\n", quiet=TRUE)
  on.exit(close(con), TRUE)

  if (!require("XML", quietly = TRUE))
    stop("read_epp_t0() requires the package 'XML'. Please install it.", call. = FALSE)
      
  obj <- xmlTreeParse(epp.xml)
  r <- xmlRoot(obj)[[1]]
  eppSetChildren.idx <- which(xmlSApply(r, xmlAttrs) == "eppSetChildren")

  t0 <- list()
  for(eppSet.idx in 1:xmlSize(r[[eppSetChildren.idx]])){

    eppSet <- r[[eppSetChildren.idx]][[eppSet.idx]][[1]]
    eppName <- xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "name")]][["string"]])
    t0[[eppName]] <- as.integer(xmlToList(eppSet[[which(xmlSApply(eppSet, xmlAttrs) == "priorT0vr")]][[1]]))
  }

  return(unlist(t0))
}


## Read subpopulation size input file
#'
#' @param filepath file path to .subp file
read_subp_file <- function(filepath){

  dat <- readLines(filepath)

  AG <- as.integer(sub("AGEGROUPS=([0-9]*)", "\\1", grep("AGEGROUPS=([0-9]*)", dat, value=TRUE)))
  startyear <- as.integer(sub("STARTYEAR=([0-9]*)", "\\1", grep("STARTYEAR=([0-9]*)", dat, value=TRUE)))
  years <- as.integer(sub("YEARS=([0-9]*)", "\\1", grep("YEARS=([0-9]*)", dat, value=TRUE)))

  startidx <- grep("DATASTART", dat)
  endidx <- grep("DATAEND", dat)

  datidx <- seq(startidx+1, endidx-1, AG+1)

  header <- setNames(read.csv(text=dat[datidx], header=FALSE),
                     c("country_code", "country", "region", "sex"))
  header$sex <- factor(header$sex, c("Male", "Female"))

  data <- lapply(datidx, function(idx) read.csv(text=dat[idx+1:AG], header=FALSE))
  data <- lapply(tapply(setNames(data, header$sex), header$region, abind::abind, along=0),
                 aperm, c(2,1,3))
  data <- lapply(data, function(x) x[,c("Male", "Female"),])
  data <- lapply(data, "dimnames<-", list(Age=0:(AG-1), Sex=c("Male", "Female"), Year=startyear+1:years-1L))

  return(data)
}


#' Read CSAVR input data
#'
#' @param pjnz file path to Spectrum PJNZ file.
read_csavr_data <- function(pjnz){

  dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  exists_dptag <- function(tag, tagcol=1){tag %in% dp[,tagcol]}
  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  yr_start <- as.integer(dpsub("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub("<FinalYear MV2>",2,4))
  proj_years <- yr_start:yr_end


  if(exists_dptag("<FitIncidenceEditorValues MV2>")){
    val <- data.frame(year = proj_years, 
                      t(sapply(dpsub("<FitIncidenceEditorValues MV2>", 2:10, 3+seq_along(proj_years)), as.numeric)),
                      row.names=proj_years)
    names(val) <- c("year", "plhiv", "plhiv_undercount", "new_cases", "new_cases_undercount", "new_cases_lag",
                    "aids_deaths", "aids_deaths_undercount", "deaths_hivp", "deaths_hivp_undercount")
    
    attr(val, "agegroup") <-  c("All ages", "Adults 15-49", "Adults 15+")[as.integer(dpsub("<IncidenceAgeGroupIndex MV>", 2, 4))+1L]
  } else
    val <- NULL

  return(val)
}
