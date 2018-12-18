#' @import data.table

aim.dir <- paste0(root,"/WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/AIM_assumptions/")
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
  dir <- paste0('/share/hiv/epp_output/gbd19/', run.name, '/')

  ## Population
  years <- dimnames(demp$basepop)[[3]]
  pop <- fread(paste0(dir, '/population_single_age/', loc, '.csv'))
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
  surv <- fread(paste0('/share/gbd/WORK/02_mortality/03_models/5_lifetables/results/hivfree_sx/locs/', loc, '_life_tables.csv'))[,V1 := NULL]
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
  asfr <- fread(paste0(dir,'/ASFR/', loc, '.csv'))
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
  tfr <- fread(paste0(dir, '/TFR/', loc, '.csv'))
  tfr <- extend.years(tfr, years)
  tfr.list <- tfr[,value]
  names(tfr.list) <- tfr$year
  demp$tfr <- tfr.list
  
  
  ## SRB
  
  ## Migration
  mig <- fread(paste0(dir, '/migration/', loc, '.csv'))
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
  years <- epp.subp[[1]]$year
  dir <- paste0('/share/hiv/epp_output/gbd19/', run.name, '/')
  # add in missing years in the future
  in.pop <- fread(paste0(dir, '/population/', loc, '.csv'))
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
  mig.data <- fread(paste0(dir, '/migration/', loc, '.csv'))
  max.mig <- copy(mig.data[year == max(year)])
  missing.dt <- rbindlist(lapply(setdiff(years, unique(mig.data$year)), function(c.year) {
    dt <- copy(max.mig)
    dt[, year := c.year]
  }))
  mig.data <- rbind(mig.data, missing.dt)
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
  ## TODO: Where are prev surveys compiled?
  surv.path <- paste0("/ihme/hiv/epp_output/gbd19/", run.name, "/prev_surveys.csv")
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

sub.prev.granular <- function(dt, loc){
  age.prev.dt <- fread('/homes/tahvif/age_prev_surveys.csv')
  age.prev.dt <- age.prev.dt[iso3 == loc]
  age.prev.dt[, agegr := paste0(age_year, '-', age_year+4)]
  age.prev.dt[,sex := ifelse(sex_id == 1, 'male', 'female')]
  ## TODO: What do these mean
  age.prev.dt[,c('used', 'deff', 'deff_approx') := c(TRUE, 2, 2)]
  age.prev.dt <- age.prev.dt[,.(year, sex, agegr, n, prev, se, used, deff, deff_approx)]
  gen.pop.dict <- c("General Population", "General population", "GP", "GENERAL POPULATION", "GEN. POPL.", "General population(Low Risk)", "Remaining Pop")
  if(length(dt) == 1) {
    gen.pop.i <- 1
  } else {
    gen.pop.i <- which(names(dt) %in% gen.pop.dict)
  }
  attr(dt[[gen.pop.i]], 'eppd')$hhs <- as.data.frame(age.prev.dt)
  return(dt)
}

sub.off.art <- function(dt, loc, k) {
  # Off-ART Mortality
  mortnoart <- fread(paste0(aim.dir, "transition_parameters/HIVmort_noART/current_draws_combined/",loc,"_mortality_par_draws.csv"))
  mortnoart <- mortnoart[draw==k,]
  mortnoart[age == '45-100', age := '45+']
  mortnoart[age == '15-25', age := '15-24']
  mortnoart[age == '25-35', age := '25-34']
  mortnoart[age == '35-45', age := '35-44']
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
  ## ???
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(7, length(dimnames(attr(dt[[1]], 'specfp')$cd4_mort)$agecat), 2))
  dimnames(replace) <- list(cd4stage = paste0(1:7), agecat = dimnames(attr(dt[[1]], 'specfp')$cd4_mort)$agecat, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(dimnames(attr(dt[[1]], 'specfp')$cd4_mort)$agecat)){
        for(c.cd4 in paste0(1:7)){
          replace[c.cd4, c.ageindex, c.sex] = as.numeric(mortnoart[age == age.index[[c.ageindex]] & cat == c.cd4, risk])
      }
    }
  }  
  
  for (n in names(dt)) {
    for(i in 1:7){
      attr(dt[[n]], 'specfp')$cd4_mort <- replace
    }
  }
  return(dt)
}

sub.on.art <- function(dt, loc, k) {
  mortart <- fread(paste0(aim.dir,"transition_parameters/HIVmort_onART_regions/DisMod/", loc,"_HIVonART.csv"))
  mortart <- melt(mortart, 
                  id = c("durationart", "cd4_category", "age", "sex","cd4_lower",
                         "cd4_upper"))

  setnames(mortart, c("variable","value","cd4_category"),c("draw","mort","cd4"))
  # mortart <- mortart[age!="55-100",]

  mortart <- mortart[draw==paste0('mort',k),]
  mortart[,sex := as.character(sex)]
  mortart[, sex := ifelse(sex == 1, 'Male', 'Female')]
  ## TODO: Find a better way to align 45-55 and 55+
  mortart[age == '45-55', age := '45+']
  mortart[age == '15-25', age := '15-24']
  mortart[age == '25-35', age := '25-34']
  mortart[age == '35-45', age := '35-44']
  mortart <- mortart[age != '55-100']
  # mortart[,age:= as.integer(sapply(strsplit(mortart[,age],'-'), function(x) {x[1]}))]
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
  mortart <- mortart[,.(durationart, age, sex, mort, cat)]
  mortart[,risk:=-1*log(1-mort)]
  mortart[, mort := NULL]
  mortart[durationart == '6to12Mo', artdur := 'ART6MOS']
  mortart[durationart == 'GT12Mo', artdur := 'ART1YR']
  mortart[durationart == 'LT6Mo', artdur := 'ART0MOS']
  mortart[, durationart := NULL]
  setnames(mortart, c('age', 'cat'), c('agecat', 'cd4stage'))
  
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(3, 7, length(dimnames(attr(dt[[1]], 'specfp')$art_mort)$agecat), 2))
  dimnames(replace) <- list(artdur = c('ART0MOS', 'ART6MOS', 'ART1YR'), cd4stage = paste0(1:7), agecat = dimnames(attr(dt[[1]], 'specfp')$art_mort)$agecat, sex = c('Male', 'Female'))
  
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(dimnames(attr(dt[[1]], 'specfp')$art_mort)$agecat)){
      for(c.dur in c('ART0MOS', 'ART6MOS', 'ART1YR')){
        for(c.cd4 in paste0(1:7)){
            replace[c.dur, c.cd4, c.ageindex, c.sex] = as.numeric(mortart[sex == c.sex & artdur == c.dur & agecat == age.index[[c.ageindex]] & cd4stage == c.cd4, risk])
        }
      }
    }
  }
  for(n in names(dt)){
    attr(dt[[n]], 'specfp')$art_mort <- replace
  }
  
  return(dt)
}

sub.cd4.prog <- function(dt, loc, k){
  progdata <- fread(paste0(aim.dir, "transition_parameters/DurationCD4cats/current_draws_combined/", loc, "_progression_par_draws.csv"))
  progdata <- progdata[order(age,cd4,draw)]
  progdata_read <- progdata[,c("age","cd4","draw","prog"), with=F]
  progdata_read <- progdata_read[,lambda:=1/prog]
  
  progdata <- progdata_read[draw==k,]
  progdata[,risk:=-1*log(1-prog)/0.1]
  progdata[,prob:=1-exp(-1*risk)]
  
  progdata[age == '45-100', age := '45+']
  progdata[age == '15-25', age := '15-24']
  progdata[age == '25-35', age := '25-34']
  progdata[age == '35-45', age := '35-44']
  progdata <- progdata[cd4=="GT500CD4", cat := 1]
  progdata <- progdata[cd4=="350to500CD4", cat := 2]
  progdata <- progdata[cd4=="250to349CD4", cat := 3]
  progdata <- progdata[cd4=="200to249CD4", cat := 4]
  progdata <- progdata[cd4=="100to199CD4", cat := 5]
  progdata <- progdata[cd4=="50to99CD4", cat := 6] 
  progdata[,risk:=-1*log(1-prob)]
  
  age.index <- list(a = '15-24', b = '15-24', c = '15-24', d = '25-34', e = '25-34', f = '35-44', g = '35-44', h = '45+', i = '45+')
  replace <- array(0, c(6, length(dimnames(attr(dt[[1]], 'specfp')$cd4_prog)$agecat), 2))
  dimnames(replace) <- list(cd4stage = paste0(1:6), agecat = dimnames(attr(dt[[1]], 'specfp')$cd4_prog)$agecat, sex = c('Male', 'Female'))
  for(c.sex in c('Male', 'Female')){
    for(c.ageindex in 1:length(dimnames(attr(dt[[1]], 'specfp')$cd4_prog)$agecat)){
      for(c.cd4 in paste0(1:6)){
        replace[c.cd4, c.ageindex, c.sex] = as.numeric(progdata[age == age.index[[c.ageindex]] & cat == c.cd4, risk])
      }
    }
  }  
  for (n in names(dt)) {
    attr(dt[[n]], 'specfp')$cd4_prog <- replace
  }	
  return(dt)
}
