
###################################################
####  function to read HIV projection outputs  ####
###################################################

read.hivproj.output <- function(specdp.file, single.age=TRUE){

  ## read .DP file
  dp <- read.csv(specdp.file, as.is=TRUE)

  dp.vers <- dp[2,1] # <General 3>: 2013, 2014 Spectrum files; <General5>: 2015 Spectrum files

  version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
  validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
  validversion <- as.numeric(dp[which(dp[,1] == "<ValidVers>")+2,4])

  ## projection parameters
  yr.start <- as.integer(dp[which(dp[,2] == "First year")+1,4])
  yr.end <- as.integer(dp[which(dp[,2] == "Final year")+1,4])
  proj.years <- yr.start:yr.end
  timedat.idx <- 4+1:length(proj.years)-1

  agegr.lab <- c(paste(0:15*5, 1:16*5, sep="-"), "80+")


  ## Number HIV+
  hivnum.tidx <- which(dp[,2] == "HIV")
  hivnum.age <- sapply(dp[hivnum.tidx+4:54, timedat.idx], as.numeric)
  rownames(hivnum.age) <- dp[hivnum.tidx+4:54, 3]
  hivnum.m <- hivnum.age[grep("Sex=1", rownames(hivnum.age)),]
  hivnum.f <- hivnum.age[grep("Sex=2", rownames(hivnum.age)),]
  dimnames(hivnum.m) <- dimnames(hivnum.f) <- list(agegr.lab, proj.years)

  ## Number new infections
  newinf.tidx <- which(dp[,2] == "New Infections")
  newinf.age <- sapply(dp[newinf.tidx+4:54, timedat.idx], as.numeric)
  rownames(newinf.age) <- dp[newinf.tidx+4:54, 3]
  newinf.m <- newinf.age[grep("Sex=1", rownames(newinf.age)),]
  newinf.f <- newinf.age[grep("Sex=2", rownames(newinf.age)),]
  dimnames(newinf.m) <- dimnames(newinf.f) <- list(agegr.lab, proj.years)

  ## Number of infant (age 0) new infections
  
  ## Total population size
  if(dp.vers == "<General 3>"){
    totpop.tidx <- which(dp[,2] == "Total Population")
    totpop.m <- sapply(dp[totpop.tidx + 1:17*7 + 6, timedat.idx], as.numeric)
    totpop.f <- sapply(dp[totpop.tidx + 1:17*7 + 8, timedat.idx], as.numeric)
    dimnames(totpop.m) <- dimnames(totpop.f) <- list(agegr.lab, proj.years)
  } else if(dp.vers == "<General5>"){
    totpop.tidx <- which(dp[,1] == "<BigPop3>")
    totpop.m <- sapply(lapply(dp[totpop.tidx + 1:81 + 1, timedat.idx], as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    totpop.f <- sapply(lapply(dp[totpop.tidx + 1:81 + 82, timedat.idx], as.numeric), tapply, c(rep(1:16, each=5), 17), sum)
    dimnames(totpop.m) <- dimnames(totpop.f) <- list(agegr.lab, proj.years)
  }

  ## ART need
  artneed.tidx <- which(dp[,2] == "Need FL")
  artneed.m <- sapply(dp[artneed.tidx+0:16*3+5, timedat.idx], as.numeric)
  artneed.f <- sapply(dp[artneed.tidx+0:16*3+6, timedat.idx], as.numeric)
  dimnames(artneed.m) <- dimnames(artneed.f) <- list(agegr.lab, proj.years)

  ## On ART
  artnum.tidx <- which(dp[,2] == "On FL")
  artnum.m <- sapply(dp[artnum.tidx+0:16*3+5, timedat.idx], as.numeric)
  artnum.f <- sapply(dp[artnum.tidx+0:16*3+6, timedat.idx], as.numeric)
  dimnames(artnum.m) <- dimnames(artnum.f) <- list(agegr.lab, proj.years)

  ## number non-aids deaths
  natdeathm.tidx <- which(dp[,2] == "Deaths - Male")
  natdeathf.tidx <- which(dp[,2] == "Deaths - Female")
  natdeaths.m <- sapply(dp[natdeathm.tidx+2:18, timedat.idx], as.numeric)
  natdeaths.f <- sapply(dp[natdeathf.tidx+2:18, timedat.idx], as.numeric)
  dimnames(natdeaths.m) <- dimnames(natdeaths.f) <- list(agegr.lab, proj.years)

  ## number AIDS deaths
  aidsdm.tidx <- which(dp[,2] == "AIDS deaths - Male")
  aidsdf.tidx <- which(dp[,2] == "AIDS deaths - Female")
  aidsdeaths.m <- sapply(dp[aidsdm.tidx+2:18, timedat.idx], as.numeric)
  aidsdeaths.f <- sapply(dp[aidsdf.tidx+2:18, timedat.idx], as.numeric)
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

  
  class(specres) <- "specres"

  return(specres)
}


######################################################
####  function to read HIV projection parameters  ####
######################################################

read.hivproj.param <- function(specdp.file){

  ## read .DP file
  dp <- read.csv(specdp.file, as.is=TRUE)

  dp.vers <- dp[2,1] # <General 3>: 2013, 2014 Spectrum files; <General5>: 2015 Spectrum files

  version <- as.numeric(dp[which(dp[,2] == "Version")+1,4])
  validdate <- dp[which(dp[,1] == "<ValidDate>")+2,3]
  validversion <- as.numeric(dp[which(dp[,1] == "<ValidVers>")+2,4])

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

  nathist.tidx <- which(dp[,1] == "<AdultTransParam2>")
  cd4initdist.tidx <- which(dp[,1] == "<DistNewInfectionsCD4>")
  infectreduc.tidx <- which(dp[,1] == "<InfectReduc>")
  adult.artnumperc.tidx <- which(dp[,1] == "<HAARTBySexPerNum>")
  adult.art.tidx <- which(dp[,1] == "<HAARTBySex>")
  adult.arteligthresh.tidx <- which(dp[,1] == "<CD4ThreshHoldAdults>")
  specpopelig.tidx <- which(dp[,1] == "<PopsEligTreat1>")
  

  ## state space dimensions
  NG <- 2
  AG <- 17
  DS <- 8
  TS <- 4

  ## projection parameters
  yr.start <- as.integer(dp[which(dp[,2] == "First year")+1,4])
  yr.end <- as.integer(dp[which(dp[,2] == "Final year")+1,4])
  proj.years <- yr.start:yr.end
  t0 <- as.numeric(dp[epidemfirstyr.tidx+2,4])
  timedat.idx <- 4+1:length(proj.years)-1

  ## scalar paramters
  relinfect.ART <- 1.0 - as.numeric(dp[infectreduc.tidx+1, 4])
  vert.trans <- 0  ## !! Spectrum vertical transmission not yet implemented

  if(dp.vers == "<General 3>"){
    fert.rat <- as.numeric(dp[which(dp[,1] == "<AIDS5>")+185, 4+0:6])
    fert.rat <- array(rep(fert.rat, length(proj.years)), c(7, length(proj.years)))
    dimnames(fert.rat) <- list(seq(15, 45, 5), proj.years)
  } else if(dp.vers == "<General5>") {
    fert.rat <- sapply(dp[hivtfr.tidx+2:8, 3+seq_along(proj.years)], as.numeric)
    dimnames(fert.rat) <- list(seq(15, 45, 5), proj.years)
  }

  ## sex/age-specific incidence ratios (time varying)
  if(dp.vers == "<General 3>"){
    inc.sexrat <- setNames(as.numeric(dp[aids5.tidx+181,timedat.idx]), proj.years) # !!! Not sure what aids5.tidx+183 (Ratio of female to male prevalence)
    inc.agerat <- array(NA, c(AG, NG, length(proj.years)), list(0:(AG-1)*5, c("Male", "Female"), proj.years))
    inc.agerat[,"Male",] <- sapply(dp[aids5.tidx+281:297,timedat.idx], as.numeric)
    inc.agerat[,"Female",] <- sapply(dp[aids5.tidx+299:315,timedat.idx], as.numeric)
  } else if(dp.vers == "<General5>"){
    inc.sexrat <- setNames(as.numeric(dp[hivsexrat.tidx+2, timedat.idx]), proj.years)
    inc.agerat <- array(NA, c(AG, NG, length(proj.years)), list(0:(AG-1)*5, c("Male", "Female"), proj.years))
    inc.agerat[,"Male",] <- sapply(dp[hivagedist.tidx+3:19,timedat.idx], as.numeric)
    inc.agerat[,"Female",] <- sapply(dp[hivagedist.tidx+21:37,timedat.idx], as.numeric)
  }

  ## hiv natural history
  cd4.initdist <- array(NA, c(DS-1, AG, NG), list(2:DS, 0:(AG-1)*5, c("Male", "Female")))
  cd4.initdist[,,"Male"] <- array(as.numeric(dp[cd4initdist.tidx+2, 4:31])/100, c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.initdist[,,"Female"] <- array(as.numeric(dp[cd4initdist.tidx+3, 4:31])/100, c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]

  cd4.prog <- array(NA, c(DS-2, AG, NG), list(2:(DS-1), 0:(AG-1)*5, c("Male", "Female")))
  cd4.prog[,,"Male"] <- array(1/as.numeric(dp[nathist.tidx+4, 4:27]), c(DS-2, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.prog[,,"Female"] <- array(1/as.numeric(dp[nathist.tidx+5, 4:27]), c(DS-2, 4))[,rep(1:4, c(5,2,2,8))]

  cd4.art.mort <- array(NA, c(TS, DS-1, AG, NG), list(1:TS, 2:DS, 0:(AG-1)*5, c("Male", "Female")))
  cd4.art.mort[1,,,"Male"] <- array(as.numeric(dp[nathist.tidx+7, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[1,,,"Female"] <- array(as.numeric(dp[nathist.tidx+8, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[2,,,"Male"] <- array(as.numeric(dp[nathist.tidx+10, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[2,,,"Female"] <- array(as.numeric(dp[nathist.tidx+11, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[3,,,"Male"] <- array(as.numeric(dp[nathist.tidx+13, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[3,,,"Female"] <- array(as.numeric(dp[nathist.tidx+14, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[4,,,"Male"] <- array(as.numeric(dp[nathist.tidx+16, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]
  cd4.art.mort[4,,,"Female"] <- array(as.numeric(dp[nathist.tidx+17, 4:31]), c(DS-1, 4))[,rep(1:4, c(5,2,2,8))]

  ## program parameters
  artnumperc.15plus <- sapply(dp[adult.artnumperc.tidx+3:4, timedat.idx], as.numeric)
  dimnames(artnumperc.15plus) <- list(c("Male", "Female"), proj.years)

  artnum.15plus <- sapply(dp[adult.art.tidx+3:4, timedat.idx], as.numeric)
  dimnames(artnum.15plus) <- list(c("Male", "Female"), proj.years)

  arteligthresh.15plus <- setNames(as.numeric(dp[adult.arteligthresh.tidx+2, timedat.idx]), proj.years)

  artelig.specpop <- setNames(dp[specpopelig.tidx+1:7,2:6], c("description", "pop", "elig", "percent", "year"))
  artelig.specpop$pop <- c("PW", "TBHIV", "DC", "FSW", "MSM", "IDU", "OTHER")
  artelig.specpop$elig <- as.logical(as.integer(artelig.specpop$elig))
  artelig.specpop$percent <- as.numeric(artelig.specpop$percent)/100
  artelig.specpop$year <- as.integer(artelig.specpop$year)
  artelig.specpop$idx <- match(as.integer(artelig.specpop$year), proj.years)
  rownames(artelig.specpop) <- artelig.specpop$pop

  ## permute indices to match model state space (year, NG, AG, DS, TS)
  fert.rat <- aperm(fert.rat, 2:1)
  inc.agerat <- aperm(inc.agerat, 3:1)
  cd4.initdist <- aperm(cd4.initdist, 3:1)
  cd4.prog <- aperm(cd4.prog, 3:1)
  cd4.art.mort <- aperm(cd4.art.mort, 4:1)
  artnum.15plus <- aperm(artnum.15plus, 2:1)
  artnumperc.15plus <- aperm(artnumperc.15plus, 2:1)

  projp <- list("yr.start"=yr.start, "yr.end"=yr.end, "t0"=t0,
                "relinfect.ART"=relinfect.ART, "vert.trans"=vert.trans,
                "fert.rat"=fert.rat, "inc.sexrat"=inc.sexrat, "inc.agerat"=inc.agerat,
                "cd4.initdist"=cd4.initdist, "cd4.prog"=cd4.prog, "cd4.art.mort"=cd4.art.mort,
                "artnumperc.15plus"=artnumperc.15plus, "artnum.15plus"=artnum.15plus,
                "arteligthresh.15plus"=arteligthresh.15plus, "artelig.specpop"=artelig.specpop)
  class(projp) <- "projp"
  attr(projp, "version") <- version
  attr(projp, "validdate") <- validdate
  attr(projp, "validversion") <- validversion

  return(projp)
}


###################################################################
####  function to read UN Population Division projection file  ####
###################################################################


read.demog.param <- function(upd.file, age.intervals = 1){

  DEM.INIT.YEAR <- 1970

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

  bp.tidx <- which(upd[,1] == "<basepop>")
  lt.tidx <- which(upd[,1] == "<lfts>")
  pasfrs.tidx <- which(upd[,1] == "<pasfrs>")
  migration.tidx <- which(upd[,1] == "<migration>")
  tfr.tidx <- which(upd[,1] == "<tfr>")
  srb.tidx <- which(upd[,1] == "<srb>")

  bp <- setNames(data.frame(upd[bp.tidx+2:1459,1:4]), upd[bp.tidx+1,1:4])
  lt <- setNames(data.frame(upd[lt.tidx+2:13121,1:12]), upd[lt.tidx+1,1:12])
  pasfrs <- setNames(data.frame(upd[pasfrs.tidx+1+1:(80*35),1:3]), upd[pasfrs.tidx+1,1:3])
  migration <- setNames(data.frame(upd[migration.tidx+1+1:(80*2*81),1:4]), upd[migration.tidx+1,1:4])
  tfr <- setNames(as.numeric(upd[tfr.tidx+1+1:80,2]), upd[tfr.tidx+1+1:80,1])
  srb <- setNames(as.numeric(upd[srb.tidx+1+1:80,2]), upd[srb.tidx+1+1:80,1])


  ## Aggregate into specified age groups

  ## population size
  basepop <- array(as.numeric(bp$value), c(81, 2, 9))
  dimnames(basepop) <- list(0:80, c("Male", "Female"), seq(1970, 2010, 5))
  basepop <- apply(basepop, 2:3, tapply, age.groups, sum)

  ## mx
  Sx <- as.numeric(lt$Sx[-(1:(2*80)*82-1)]) # 80+ age group given twice
  dim(Sx) <- c(81, 2, 80)
  dimnames(Sx) <- list(0:80, c("Male", "Female"), 1970:2049)
  Sx <- apply(Sx, 2:3, tapply, age.groups, prod)
  mx <- -sweep(log(Sx), 1, age.intervals, "/")

  ## asfr
  asfd <- array(as.numeric(pasfrs$value), c(35, 80))
  dimnames(asfd) <- list(15:49, 1970:2049)
  asfr <- sweep(asfd, 2, tfr, "*")
  asfr <- apply(asfr, 2, tapply, age.groups[16:50], mean)

  asfd <- apply(asfd, 2, tapply, age.groups[16:50], sum)

  ## migration
  netmigr <- array(as.numeric(migration$value), c(81, 2, 80))
  netmigr <- apply(netmigr, 2:3, tapply, age.groups, sum)

  ## ## permute indices to match model state space (year, NG, AG)
  ## basepop <- aperm(basepop, 3:1)
  ## mx <- aperm(mx, 3:1)
  ## Sx <- aperm(Sx, 3:1)
  ## asfr <- aperm(asfr, 2:1)
  ## asfd <- aperm(asfd, 2:1)
  ## netmigr <- aperm(netmigr, 3:1)

  demp <- list("basepop"=basepop, "mx"=mx, "Sx"=Sx, "asfr"=asfr, "tfr"=tfr, "asfd"=asfd, "srb"=srb, "netmigr"=netmigr)
  class(demp) <- "demp"

  return(demp)
}
