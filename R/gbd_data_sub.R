#' @import data.table
extend.years <- function(dt, years){
  dt <- as.data.table(dt)
  if('year_id' %in% names(dt)){setnames(dt, 'year_id', 'year')}
  dt <- dt[as.integer(year) %in% as.integer(years)]
  while(max(dt$year) < max(as.integer(years))){
    dt.ext <- dt[year == max(dt$year)]
    dt.ext[, year := year + 1]
    dt <- rbind(dt, dt.ext, use.names = T)
  }  
  return(dt)
}

sub.pop.params.demp <- function(demp, loc, k){
  ## Population
  years <- dimnames(demp$basepop)[[3]]
  pop <- fread(paste0('/Users/tahvif/Documents/code/eppasm/gbd-data/', loc, '_single_age_pop.csv'))
  pop <- extend.years(pop, years)
  pop[,age := ifelse(age_group_id == 28, 0, age_group_id - 48)]
  pop <- pop[age %in% 0:80]
  pop[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  pop <- dcast.data.table(pop[,.(age, sex, year, population)], age + year ~ sex, value.var = 'population')
  for(i in 1:length(years)){
    pop.year <- as.matrix(pop[year == as.integer(years[i]),.(Male, Female)])
    rownames(pop.year) <- 0:80
    demp$basepop[,,i] <- pop.year
  }
  
  ## Survival
  surv <- fread(paste0('/Users/tahvif/Documents/code/eppasm/gbd-data/', loc, '_life_tables.csv'))[,V1 := NULL]
  surv <- melt(surv, id.vars = c('sex', 'year', 'age'))
  surv[, variable := as.integer(gsub('px', '', variable))]
  surv <- surv[variable == k]
  surv[,variable := NULL]
  surv[sex == 'male', sex := 'Male']
  surv[sex == 'female', sex := 'Female']
  surv <- dcast.data.table(surv, year + age ~ sex)
  surv = extend.years(surv, years)
  surv <- surv[age %in% 0:80]
  for(i in 1:length(years)){
    sx.year = as.matrix(surv[year == as.integer(years[i]), .(Male, Female)])
    rownames(sx.year) <- 0:80
    demp$Sx[,,i] <- sx.year
  }
  
  ## mx
  mx <- surv
  mx[, Female := -log(Female)]
  mx[, Male := -log(Male)]
  for(i in 1:length(years)){
    mx.year = as.matrix(mx[year == as.integer(years[i]), .(Male, Female)])
    rownames(mx.year) <- 0:80
    demp$mx[,,i] <- sx.year
  }
  
  ## ASFR
  asfr <- fread('/Users/tahvif/Documents/code/eppasm/gbd-data/MWI_ASFR.csv')
  asfr <- extend.years(asfr, years)
  ## Copy 5-year asfr
  for(c.age in 15:49){
    if(!c.age %in% asfr$age){
      asfr.ext <- asfr[age == (c.age - c.age%%5)]
      asfr.ext[, age := c.age]
      asfr <- rbind(asfr, asfr.ext)
    }
  }
  asfr <- as.matrix(dcast.data.table(asfr, age~year))
  asfr <- asfr[,2:(length(years) + 1)]
  rownames(asfr) <- 15:49
  demp$asfr <- asfr
  
  ## TFR
  tfr <- fread(paste0('/Users/tahvif/Documents/code/eppasm/gbd-data/', loc, '_TFR.csv'))
  tfr <- extend.years(tfr, years)
  tfr.list <- tfr[,value]
  names(tfr.list) <- tfr$year
  demp$tfr <- tfr.list
  
  
  ## SRB
  
  ## Migration
  mig <- fread(paste0('/Users/tahvif/Documents/code/eppasm/gbd-data/', loc, '_migration.csv'))
  setnames(mig, 'sex', 'sex_id')
  mig[, sex := ifelse(sex_id == 1, 'Male', 'Female')]
  mig <- dcast.data.table(mig[,.(year, age, value, sex)], year + age ~ sex)
  mig <- extend.years(mig, years)
  for(i in 1:length(years)){
    mig.year = as.matrix(mig[year == as.integer(years[i]), .(Male, Female)])
    rownames(mig.year) <- 0:80
    demp$netmigr[,,i] <- mig.year
  }
  return(demp)
}

sub.pop.params.epp <- function(epp.subp, epp.input, loc) {
  ## Load central functions
  # loc.id <- loc.table[ihme_loc_id==loc, location_id]
  years <- epp.subp[[1]]$year
  path <- paste0("/Users/tahvif/Documents/code/eppasm/gbd-data/", loc, "_pop_all_age.csv")
  if(file.exists(path)) {
    in.pop <- fread(path)
    in.pop <- in.pop[year_id %in% years]
  } else {
    if(!"get_population" %in% ls()) {
      source(paste0(root, "temp/central_comp/libraries/current/r/get_population.R"))
    }	
    in.pop <- get_population(age_group_id = 8:14, location_id = loc.id, year_id = years, sex_id = 1:2, location_set_id = 21)
  }
  
  # add in missing years in the future
  max.pop <- copy(in.pop[year_id == max(year_id)])
  missing.dt <- rbindlist(lapply(setdiff(years, unique(in.pop$year_id)), function(year) {
    dt <- copy(max.pop)
    dt[, year_id := year]
  }))
  bound.pop <- rbind(in.pop, missing.dt)
  bound.pop <- bound.pop[age_group_id %in% 8:14]
  both.pop <- bound.pop[, .(population = sum(population)), by = .(age_group_id, year_id)]
  pop15 <- both.pop[age_group_id == 8, population] / 5
  pop50 <- both.pop[age_group_id == 14, population] / 5
  pop15to49 <- both.pop[, .(population = sum(population)), by = .(year_id)]$population
  
  ## TODO: Is migration 15-49 sum?
  mig.data <- fread('/Users/tahvif/Documents/code/eppasm/gbd-data/MWI_migration.csv')
  mig.data <- mig.data[age %in% 15:49 & year %in% both.pop$year_id]
  mig.data <- mig.data[,.(value = sum(value)), by = .(year)]$value
  
  for (pop in names(epp.subp)) {
    temp.dt <- epp.subp[[pop]]
    if(pop=="subpops") {
      for (subpop in names(temp.dt)){       
        epp.subp[[pop]][[subpop]]$pop15to49 <- pop15to49
        epp.subp[[pop]][[subpop]]$pop15 <- pop15
        epp.subp[[pop]][[subpop]]$pop50 <- pop50
        epp.subp[[pop]][[subpop]]$netmigr <- mig.data
      } 
    } else {
      epp.subp[[pop]]$pop15to49 <- pop15to49
      epp.subp[[pop]]$pop15 <- pop15
      epp.subp[[pop]]$pop50 <- pop50
      epp.subp[[pop]]$netmigr <- mig.data
    }
  }
  epp.input[['epp.pop']]$pop15to49 <- pop15to49
  epp.input[['epp.pop']]$pop15 <- pop15
  epp.input[['epp.pop']]$pop50 <- pop50
  epp.input[['epp.pop']]$netmigr <- mig.data
  return(list(epp.subp = epp.subp, epp.input = epp.input))
}

sub.prev <- function(loc, dt){
  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
  if(length(dt) == 1) {
    gen.pop.i <- 1
  } else {
    gen.pop.i <- which(names(dt) %in% gen.pop.dict)
  }
  surv.path <- paste0("/Users/tahvif/Documents/code/eppasm/gbd-data/prev_surveys.csv")
  data4 <- fread(surv.path)[iso3 == loc]
  data4[,c("iso3", "int_year", "nid") := NULL]
  
  if(nrow(data4) > 0) {
    data4[,used:=TRUE]
    data4[prev==0,used:=FALSE]
    data4[,W.hhs:=qnorm(prev)]
    data4[,v.hhs:=2*pi*exp(W.hhs^2)*se^2]
    data4[,sd.W.hhs := sqrt(v.hhs)]
    ## TODO: Correct index?
    data4[,idx := year - (attr(dt[[loc]], 'eppfp')$tsEpidemicStart-1.5)]
    data4[, agegr := '15-49']
    data4[, sex := 'both']
    ## TODO: what is deff
    data4[, deff := 2]
    data4[, deff_approx := 2]
    if(grepl("ZAF", loc)) {
      drop.years <- unique(c(data4$year,data4$year+1,data4$year-1))
      data4 <- rbind(data4,attr(dt[[1]], 'likdat')$hhslik.dat[! attr(dt[[1]], 'likdat')$hhslik.dat$year %in% drop.years, ])
    } 
    data4 <- data4[order(data4$year),]
    
    attr(dt[[gen.pop.i]], 'eppd')$hhs <- as.data.frame(data4[, .(year, sex, agegr, n, prev, se, used, deff, deff_approx)])
    if(all(!attr(dt[[gen.pop.i]], "eppd")$anc.used)) {
      attr(dt[[gen.pop.i]], 'likdat')$hhslik.dat <- data4
    } else {
      attr(dt[[gen.pop.i]], 'likdat') <- prepare_likdat(attr(dt[[gen.pop.i]], 'eppd'), attr(dt[[gen.pop.i]], "specfp"))      
    }
  } else { 
    print(paste0("No surveys for ",loc))
  }
  return(dt)
}

sub.off.art <- function(dt, loc, k) {
  # Off-ART Mortality
  mortnoart <- fread(paste0(aim.dir, "transition_parameters/HIVmort_noART/current_draws_combined/",loc,"_mortality_par_draws.csv"))
  mortnoart[,draw:=rank(-mort,ties.method="first"),by=c("age","cd4")]
  mortnoart <- mortnoart[order(age,cd4,draw)]
  mortnoart_read <- mortnoart[,c("age","cd4","draw","mort"), with=F]
  mortnoart <- mortnoart_read[draw==k,]
  mortnoart[,age:= as.integer(sapply(strsplit(mortnoart[,age],'-'), function(x) {x[1]}))]
  mortnoart[,risk:=-1*log(1-mort)/0.1]
  mortnoart[,prob:=1-exp(-1*risk)]
  
  cd4_cats <- unique(mortnoart[,cd4])
  cd4_vars <- data.table(cd4=cd4_cats)
  
  mortnoart <- mortnoart[cd4=="GT500CD4", cat := 1]
  mortnoart <- mortnoart[cd4=="350to500CD4", cat := 2]
  mortnoart <- mortnoart[cd4=="250to349CD4", cat := 3]
  mortnoart <- mortnoart[cd4=="200to249CD4", cat := 4]
  mortnoart <- mortnoart[cd4=="100to199CD4", cat := 5]
  mortnoart <- mortnoart[cd4=="50to99CD4", cat := 6] 
  mortnoart <- mortnoart[cd4=="LT50CD4", cat := 7] 
  mortnoart[,risk:=-1*log(1-prob)]
  mortnoart <- mortnoart[,.(age,risk,cat)]
  ## TODO: Why do these seem so high compared to pjnz?
  mortnoart <- dcast.data.table(mortnoart, cat~age, value.var = 'risk')
  mortnoart <- data.frame(mortnoart)
  names(mortnoart) <- NULL
  mugbd <- as.matrix(mortnoart)[,2:5]
  
  for (n in names(dt)) {
    for(i in 1:7){
      attr(dt[[n]], 'eppfp')$cd4artmort[i,1:4] <- mugbd[i, 1:4]
    }
  }
  return(dt)
}

## TODO
sub.on.art <- function(dt, loc, k) {
  mortart <- fread(paste0("/Users/tahvif/Documents/code/eppasm/gbd-data/MWI_HIVonART.csv"))
  mortart <- melt(mortart, 
                  id = c("durationart", "cd4_category", "age", "sex","cd4_lower",
                         "cd4_upper"))

  setnames(mortart, c("variable","value","cd4_category"),c("draw","mort","cd4"))
  mortart <- mortart[age!="55-100",]

  mortart <- mortart_read[draw==k,]
  mortart[,age:= as.integer(sapply(strsplit(mortart[,age],'-'), function(x) {x[1]}))]
  mortart[,sex:=as.integer(sex)]
  cd4_cats <- unique(mortart[,cd4])
  durat_cats <- unique(mortart[,durationart])
  cd4_vars <- expand.grid(durationart=durat_cats, cd4=cd4_cats)
  mortart <- mortart[cd4=="ARTGT500CD4", cat := 1]
  mortart <- mortart[cd4=="ART350to500CD4", cat := 2]
  mortart <- mortart[cd4=="ART250to349CD4", cat := 3]
  mortart <- mortart[cd4=="ART200to249CD4", cat := 4]
  mortart <- mortart[cd4=="ART100to199CD4", cat := 5]
  mortart <- mortart[cd4=="ART50to99CD4", cat := 6] 
  mortart <- mortart[cd4=="ARTLT50CD4", cat := 7]
  mortart[,risk:=-1*log(1-mort)]
  mortart <- mortart[, setattr(as.list(risk), 'names', cat), by=c("year","durationart")]
  mortart <- mortart[order(year, durationart)]
  # mortart <- mortart[age!="55-100",]
  
  mortart1 <- mortart[durationart=="LT6Mo",]
  mortart2 <- mortart[durationart=="6to12Mo",]
  mortart3 <- mortart[durationart=="GT12Mo",]
  
  alpha1gbd <- as.matrix(data.frame(mortart1[,c("1","2","3","4","5","6", "7"), with=F]))
  alpha2gbd <- as.matrix(data.frame(mortart2[,c("1","2","3","4","5","6", "7"), with=F]))
  alpha3gbd <- as.matrix(data.frame(mortart3[,c("1","2","3","4","5","6", "7"), with=F]))
  
  alpha1 <- as.vector(t(alpha1gbd))
  alpha2 <- as.vector(t(alpha2gbd))
  alpha3 <- as.vector(t(alpha3gbd))
  for (n in names(dt)) {
    attr(dt[[n]], 'eppfp')$cd4artmort[,1] <- alpha1
    attr(dt[[n]], 'eppfp')$cd4artmort[,2] <- alpha2
    attr(dt[[n]], 'eppfp')$cd4artmort[,3] <- alpha3
    
  }
  return(dt)
}
