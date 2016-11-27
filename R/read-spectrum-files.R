
###################################################
####  function to read HIV projection outputs  ####
###################################################

read_hivproj_output <- function(pjnz, single.age=TRUE){

  ## read .DP file
  dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  dp.vers <- dp[2,1] # <General 3>: 2013, 2014 Spectrum files; <General5>: 2015 Spectrum files
  if(!dp.vers %in% c("<General 3>", "<General5>"))
    dp.vers <- "Spectrum2016"

  if(dp.vers %in% c("<General 3>", "<General5>")){
    version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
    validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
    validversion <- dp[which(dp[,1] == "<ValidVers>")+2,4]
  } else if(dp.vers == "Spectrum2016") {
    version <- as.numeric(dpsub("<VersionNum MV>", 3, 4))
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
  } else {
    totpop.m <- sapply(lapply(dpsub("<BigPop MV>", 3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f <- sapply(lapply(dpsub("<BigPop MV>", 81+3:83, timedat.idx), as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
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
  } else {
    natdeaths.m <- sapply(dpsub("<Deaths MV>", 5:21, timedat.idx), as.numeric)
    natdeaths.f <- sapply(dpsub("<Deaths MV>", 24:40, timedat.idx), as.numeric)
  }
  dimnames(natdeaths.m) <- dimnames(natdeaths.f) <- list(agegr.lab, proj.years)

  ## number AIDS deaths
  if(dp.vers %in% c("<General 3>", "<General5>")){
    aidsdeaths.m <- sapply(dpsub("AIDS deaths - Male", 2:18, timedat.idx, 2), as.numeric)
    aidsdeaths.f <- sapply(dpsub("AIDS deaths - Female", 2:18, timedat.idx, 2), as.numeric)
  } else {
    aidsdeaths.m <- sapply(dpsub("<AIDSDeaths MV>", 3:19, timedat.idx), as.numeric)
    aidsdeaths.f <- sapply(dpsub("<AIDSDeaths MV>", 22:38, timedat.idx), as.numeric)
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

  ## if(single.age){
  ##   ## !!! WORKING HERE
  ##   <HIVBySingleAge>
  ## }

  ## "<HIVBySingleAge MV>"
  ## "<DeathsByAge MV>"


  class(specres) <- "specres"

  return(specres)
}


######################################################
####  function to read HIV projection parameters  ####
######################################################

read_hivproj_param <- function(pjnz){

  ## read .DP file
  dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  dpsub <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }


  dp.vers <- dp[2,1] # <General 3>: 2013, 2014 Spectrum files; <General5>: 2015 Spectrum files
  if(!dp.vers %in% c("<General 3>", "<General5>"))
    dp.vers <- "Spectrum2016"

  if(dp.vers %in% c("<General 3>", "<General5>")){
    version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
    validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
    validversion <- as.numeric(dp[which(dp[,1] == "<ValidVers>")+2,4])
  } else if(dp.vers == "Spectrum2016") {
    version <- as.numeric(dpsub("<VersionNum MV>", 3, 4))
    validdate <- dpsub("<ValidDate MV>",2,3)
    validversion <- dpsub("<ValidVers MV>",2,4)
  }


  ## find tag indexes (tidx)
  if(dp.vers == "<General 3>"){
    aids5.tidx <- which(dp[,1] == "<AIDS5>")
    epidemfirstyr.tidx <- aids5.tidx
  } else if(dp.vers == "<General5>"){
    epidemfirstyr.tidx <- which(dp[,1] == "<FirstYearOfEpidemic>")
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
    t0 <- as.numeric(dp[epidemfirstyr.tidx+2,4])
  } else if(dp.vers == "Spectrum2016"){
    yr_start <- as.integer(dpsub("<FirstYear MV>",3,4))
    yr_end <- as.integer(dpsub("<FinalYear MV>",3,4))
    t0 <- as.numeric(dpsub("<FirstYearOfEpidemic MV>",2,4))
  }
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1
    
  
  ## scalar paramters
  if(dp.vers %in% c("<General 3>", "<General5>")){
     relinfectART <- 1.0 - as.numeric(dp[infectreduc.tidx+1, 4])
  } else if(dp.vers == "Spectrum2016"){
    relinfectART <- 1.0 - as.numeric(dpsub("<AdultInfectReduc MV>",2,4))
  }

  if(dp.vers == "<General 3>"){
    fert_rat <- as.numeric(dp[which(dp[,1] == "<AIDS5>")+185, 4+0:6])
    fert_rat <- array(rep(fert_rat, length(proj.years)), c(7, length(proj.years)))
  } else if(dp.vers == "<General5>") {
    fert_rat <- sapply(dp[hivtfr.tidx+2:8, 3+seq_along(proj.years)], as.numeric)
  } else if(dp.vers == "Spectrum2016") {
    fert_rat <- sapply(dpsub("<HIVTFR MV>", 2:8, timedat.idx), as.numeric)
  }
  dimnames(fert_rat) <- list(seq(15, 45, 5), proj.years)


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
  }


  ## hiv natural history
  cd4_initdist <- array(NA, c(DS, 4, NG), list(1:DS, c("15-24", "25-34", "35-44", "45+"), c("Male", "Female")))
  cd4_prog <- array(NA, c(DS-1, 4, NG), list(1:(DS-1), c("15-24", "25-34", "35-44", "45+"), c("Male", "Female")))
  cd4_mort <- array(NA, c(DS, 4, NG), list(1:DS, c("15-24", "25-34", "35-44", "45+"), c("Male", "Female")))
  art_mort <- array(NA, c(TS, DS, 4, NG), list(c("ART0MOS", "ART6MOS", "ART1YR"), 1:DS, c("15-24", "25-34", "35-44", "45+"), c("Male", "Female")))

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
    
  } else if(dp.vers == "Spectrum2016") {
    cd4_initdist[,,"Male"] <- array(as.numeric(dpsub("<AdultDistNewInfectionsCD4 MV>", 3, 4:31))/100, c(DS, 4))
    cd4_initdist[,,"Female"] <- array(as.numeric(dpsub("<AdultDistNewInfectionsCD4 MV>", 4, 4:31))/100, c(DS, 4))

    ## Note: CD4 progression array has DS values, but should only be DS-1. Not sure what the last one is.
    cd4_prog[,,"Male"] <- array(as.numeric(dpsub("<AdultAnnRateProgressLowerCD4 MV>", 3, 4:31)), c(DS, 4))[1:(DS-1),]
    cd4_prog[,,"Female"] <- array(as.numeric(dpsub("<AdultAnnRateProgressLowerCD4 MV>", 3, 4:31)), c(DS, 4))[1:(DS-1),]

    cd4_mort[,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4NoART MV>", 3, 4:31)), c(DS, 4))
    cd4_mort[,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4NoART MV>", 3, 4:31)), c(DS, 4))
    
    art_mort[1,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV>", 3, 4:31)), c(DS, 4))
    art_mort[1,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART0to6 MV>", 4, 4:31)), c(DS, 4))
    art_mort[2,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV>", 3, 4:31)), c(DS, 4))
    art_mort[2,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithART7to12 MV>", 4, 4:31)), c(DS, 4))
    art_mort[3,,,"Male"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV>", 3, 4:31)), c(DS, 4))
    art_mort[3,,,"Female"] <- array(as.numeric(dpsub("<AdultMortByCD4WithARTGt12 MV>", 4, 4:31)), c(DS, 4))
  }

  
  ## program parameters
  if(dp.vers %in% c("<General 3>", "<General5>")){
    art15plus_numperc <- sapply(dp[adult.artnumperc.tidx+3:4, timedat.idx], as.numeric)
    art15plus_num <- sapply(dp[adult.art.tidx+3:4, timedat.idx], as.numeric)
    art15plus_eligthresh <- setNames(as.numeric(dp[adult.arteligthresh.tidx+2, timedat.idx]), proj.years)
    artelig_specpop <- setNames(dp[specpopelig.tidx+1:7,2:6], c("description", "pop", "elig", "percent", "year"))
  } else if(dp.vers == "Spectrum2016") {
    art15plus_numperc <- sapply(dpsub("<HAARTBySexPerNum MV>", 4:5, timedat.idx), as.numeric)
    art15plus_num <- sapply(dpsub("<HAARTBySex MV>", 4:5, timedat.idx), as.numeric)
    art15plus_eligthresh <- setNames(as.numeric(dpsub("<CD4ThreshHoldAdults MV>", 2, timedat.idx)), proj.years)
    artelig_specpop <- setNames(dpsub("<PopsEligTreat MV>", 3:9, 2:6), c("description", "pop", "elig", "percent", "year"))
  }
    
  dimnames(art15plus_numperc) <- list(c("Male", "Female"), proj.years)
  dimnames(art15plus_num) <- list(c("Male", "Female"), proj.years)

  artelig_specpop$pop <- c("PW", "TBHIV", "DC", "FSW", "MSM", "IDU", "OTHER")
  artelig_specpop$elig <- as.logical(as.integer(artelig_specpop$elig))
  artelig_specpop$percent <- as.numeric(artelig_specpop$percent)/100
  artelig_specpop$year <- as.integer(artelig_specpop$year)
  artelig_specpop$idx <- match(as.integer(artelig_specpop$year), proj.years)
  rownames(artelig_specpop) <- artelig_specpop$pop

  projp <- list("yr_start"=yr_start, "yr_end"=yr_end, "t0"=t0,
                "relinfectART"=relinfectART,
                "fert_rat"=fert_rat, "incrr_sex"=incrr_sex, "incrr_age"=incrr_age,
                "cd4_initdist"=cd4_initdist, "cd4_prog"=cd4_prog, "cd4_mort"=cd4_mort, "art_mort"=art_mort,
                "art15plus_numperc"=art15plus_numperc, "art15plus_num"=art15plus_num,
                "art15plus_eligthresh"=art15plus_eligthresh, "artelig_specpop"=artelig_specpop)
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
  dimnames(basepop) <- list(unique(bp$age), c("Male", "Female"), unique(bp$year))
  basepop <- apply(basepop, 2:3, tapply, age.groups, sum)

  ## mx
  years <- unique(lt$year)
  nyears <- length(years)
  Sx <- as.numeric(lt$Sx[-(1:(2*nyears)*82-1)]) # 80+ age group given twice
  dim(Sx) <- c(81, 2, nyears)
  dimnames(Sx) <- list(0:80, c("Male", "Female"), years)
  Sx <- apply(Sx, 2:3, tapply, age.groups, prod)
  mx <- -sweep(log(Sx), 1, age.intervals, "/")

  ## asfr
  asfd <- array(as.numeric(pasfrs$value), c(35, nyears))
  dimnames(asfd) <- list(15:49, years)
  asfr <- sweep(asfd, 2, tfr, "*")
  asfr <- apply(asfr, 2, tapply, age.groups[16:50], mean)

  asfd <- apply(asfd, 2, tapply, age.groups[16:50], sum)

  ## migration
  netmigr <- array(as.numeric(migration$value), c(81, 2, nyears))
  dimnames(netmigr) <- list(0:80, c("Male", "Female"), years)
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

read_specdp_demog_param <- function(pjnz){

  dpfile <- grep(".DP$", unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- read.csv(unz(pjnz, dpfile), as.is=TRUE)

  version <- paste("Spectrum", dp[which(dp[,1] == "<ValidVers MV>")+2, 4])

  ## projection parameters
  yr_start <- as.integer(dp[which(dp[,1] == "<FirstYear MV>")+3,4])
  yr_end <- as.integer(dp[which(dp[,1] == "<FinalYear MV>")+3,4])
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1


  ## population size
  bp.tidx <- which(dp[,1] == "<BigPop MV>")
  basepop <- dp[bp.tidx+3+0:161, timedat.idx]
  basepop <- array(as.numeric(unlist(basepop)), c(81, 2, length(proj.years)))
  dimnames(basepop) <- list(0:80, c("Male", "Female"), proj.years)

  ## mx
  sx.tidx <- which(dp[,1] == "<SurvRate MV>")
  Sx <- dp[sx.tidx+3+c(0:79,81, 83+0:79, 83+81), timedat.idx]
  Sx <- array(as.numeric(unlist(Sx)), c(81, 2, length(proj.years)))
  dimnames(Sx) <- list(0:80, c("Male", "Female"), proj.years)
  mx <- -log(Sx)

  ## asfr
  tfr.tidx <- which(dp[,1] == "<TFR MV>")
  asfd.tidx <- which(dp[,1] == "<ASFR MV>")

  tfr <- setNames(as.numeric(dp[tfr.tidx + 2, timedat.idx]), proj.years)
  asfd <- sapply(dp[asfd.tidx + 3:9, timedat.idx], as.numeric)/100
  asfd <- apply(asfd / 5, 2, rep, each=5)
  dimnames(asfd) <- list(15:49, proj.years)
  asfr <- sweep(asfd, 2, tfr, "*")

  ## srb
  srb.tidx <- which(dp[,1] == "<SexBirthRatio MV>")
  srb <- setNames(as.numeric(dp[srb.tidx + 2, timedat.idx]), proj.years)

  ## migration
  migrrate.tidx <- which(dp[,1] == "<MigrRate MV>")
  migaged.tidx <- which(dp[,1] == "<MigrAgeDist MV>")

  totnetmig <- sapply(dp[migrrate.tidx+c(5,8), timedat.idx], as.numeric)

  ## note: age=0 is empty in DP file, inputs start in age=1
  netmigagedist <- sapply(dp[migaged.tidx+6+c(1:17*2, 37+1:17*2), timedat.idx], as.numeric) / 100
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
  dimnames(netmigr) <- list(0:80, c("Male", "Female"), proj.years)


  demp <- list("basepop"=basepop, "mx"=mx, "Sx"=Sx, "asfr"=asfr, "tfr"=tfr, "asfd"=asfd, "srb"=srb, "netmigr"=netmigr)
  class(demp) <- "demp"
  attr(demp, "version") <- version

  return(demp)
}
