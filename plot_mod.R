

# rm(list=ls())
# devtools::install_github("mrc-ide/eppasm@new-master-csvar")
# devtools::install_github("jeffeaton/epp@deepa-dev")
# devtools::load_all("~/Documents/eppasm")

mod = locfit_aggr
mod = locfit$`Clients des travailleurs du sexe`
pjnz <- system.file("extdata/testpjnz", "senegal_2021_05_04.pjnz", package="eppasm")


pjnz = paste0("C:/Users/DeepaJahagirdar/Avenir Health Dropbox/Avenir Shared Drive/DataSets/UNAIDS/2021 Estimates/Spectrum Final Files/WCA/Niger_2021_vf.pjnz")

entry_rates = read.csv(system.file("extdata/entry_rates.csv", package = "eppasm"))
pops = epp::read_spt(pjnz)$subpops
entrant_prev = readRDS(system.file("extdata/entrantprev.RDS", package="eppasm"))
loc <- prepare_spec_fit(pjnz, proj.end=2021.5, popupdate = FALSE, popadjust = FALSE, old=TRUE)
loc <- lapply(loc, function(x) { ##Not sure why these broke, will figure out later
  attr(x, "specfp")$entrantprev <- entrant_prev
  attr(x,"specfp")$popadjust  <- FALSE
  return(x)})

library(dplyr)
library(plyr)
library(data.table)
library(ggplot2)

age_key = data.table(age = seq(15,81,by= 1))
age_group1 = c(1,1,2,2,2,
               3,3,3,3,3,
               4,4,4,4,4,
               5,5,5,5,5,
               6,6,6,6,6,
               7,7,7,7,7,
               8,8,8,8,8,rep(9,32))
age_key[,age_group := age_group1]
age_groups = 1:66
year1 = c(1,10,20,30,40,50)


ff = alply(attr(mod,"infections")[,,year1],3)
ff = lapply(ff, function(x) prop.table(x))

ff = lapply(ff, function(x) melt(x))
for(i in  1:length(year1)){
  ff[[i]]$year <- year1[i]
}

year_dat = rbindlist(ff)
setnames(year_dat, c("Var1","Var2","value"), c("age","sex","value"))
year_dat[,year := year + 1970]
ggplot(year_dat, aes(age + 14, value, fill=year)) + 
  geom_bar(stat="identity",position = "dodge",alpha=0.95) + 
  theme_bw() + ylab(lab) + xlab("Age") + ggtitle(title)  


get_mod <- function(fit){
  
  if("resample" %in% names(fit)){
    theta_ur <- apply(fit$resample,2,mean)
  } else {
    theta_ur <- fit$par
  }
  
  fit$fp <- prepare_rhybrid(fit$fp)
  fit$fp$logitiota <- TRUE
  fit$fp$incidmod <- "eppspectrum"
  param <- fnCreateParam(theta_ur, fit$fp)
  fp_par <- update(fit$fp, list = param)
  
  mod <- simmod.specfp(fp_par, VERSION = "R")

  return(mod)
}

## Age distribution plots
plot_age_dist <- function(mod, incidence=FALSE,population, lab, title, sex_both=FALSE, plhiv = FALSE){

  if(plhiv){
    ff = alply(attr(mod,population)[,,2,year1],3)
    ff = lapply(ff, function(x) prop.table(x))
  } else {
    ff = alply(attr(mod,population)[,,,year1],4)
    ff = lapply(ff, function(x) prop.table(apply(x,1:2,sum)))
  }
  
  if(incidence){
    ff = alply(attr(mod,"infections")[,,year1],3)
    ff = lapply(ff, function(x) prop.table(x))
  }
  
  ff = lapply(ff, function(x) melt(x))
  for(i in  1:length(year1)){
    ff[[i]]$year <- year1[i]
  }
  
  year_dat = rbindlist(ff)
  setnames(year_dat, c("Var1","Var2","value"), c("age","sex","value"))
  year_dat[,year := year + 1970]
  
  if(sex_both)
    gg = ggplot(year_dat, aes(age + 14, value, fill=year)) + 
    geom_bar(stat="identity",position = "dodge",alpha=0.8) + 
    theme_bw() + ylab(lab) + xlab("Age") + facet_wrap(~sex) + ggtitle(title)
  else
    gg = ggplot(year_dat, aes(age + 14, value, fill=year)) + 
    geom_bar(stat="identity",position = "dodge",alpha=0.95) + 
    theme_bw() + ylab(lab) + xlab("Age") + ggtitle(title)  
  
  return(gg)
  
}

apply(attr(mod,"infections"))

group = "Travailleurs du sexe"
nn1 = readRDS("~/Documents/tests/SEN_test19fitmod_all.RDS")
nn2 = readRDS("~/Documents/tests/SEN_test17fitmod_all.RDS")

fit1 = nn1[[group]]
fit2 = nn2[[group]]

population = "t.allpop"
plhiv = F
incidence = T
lab = "proportion of population"
title = "Age distribution of infections - Turnover FSW"
gg1 = plot_age_dist(get_mod(fit1),incidence = T, population = population, lab = lab, title = title, plhiv=plhiv)

plhiv = F
incidence = T
title = "Age distribution of infections - FSW, higher turnover"
gg2 = plot_age_dist(get_mod(fit2),population = population, lab = lab, title = title, plhiv=plhiv)

gridExtra::grid.arrange(gg1,gg2)

##Compare to data & EPP results
devtools::load_all("~/Documents/eppasm")
group = "Pop masculine restante"
fit1 = readRDS("~/HIV/tests/2023-04-09-NIGER.RDS")[[group]]
pjnz1 = paste0("C:/Users/DeepaJahagirdar/Avenir Health Dropbox/Avenir Shared Drive/DataSets/UNAIDS/2021 Estimates/Spectrum Final Files/WCA/Niger_2021_vf.pjnz")
pop = epp::read_spt(pjnz1)
epp = data.table(pop[[group]])
epp[,model := "epp"]
epp[,year := 1:nrow(epp)]
epp = melt(epp, id.vars = c("year","model"))


get_dat = function(fit, model1){
  eppasm = data.table(prev = attr(get_mod(fit),"prev15to49"), 
                      incid = attr(get_mod(fit),"incid15to49"),
                      pop = attr(get_mod(fit),"pop"),
                      model = model1,
                      year = 1:length(attr(get_mod(fit),"pop")))
  eppasm = melt(eppasm, id.vars = c("year","model"))
  return(eppasm)
}


compare = rbind(
                get_dat(fit1,"eppasm"),
                epp, fill=TRUE,use.names=TRUE)
compare[,year := year + 1969]
compare[,variable := factor(variable, levels = c("prev","incid","pop"),
                            labels = c("Prevalence","Incidence","Population"))]
ggplot(compare[year <= 2020], aes(year, value, col = model)) + geom_line() +
  facet_wrap(~variable, scales = "free_y") + theme_bw() + ylab(NULL)
anc = data.table(attr(loc[[group]],"eppd")$ancsitedat)[,.(year,prev,type="anc",se=NA)]
hhs = data.table(attr(loc[[group]],"eppd")$hhs)[,.(year,prev,se,type="hhs")]
prev_dat = rbind(anc,hhs)



epp[,year := year + 1969]
nat_dat = data.table(year = 1970:2021,
                     value= c(attr(locfit_aggr[[1]][[1]],"incid15to49"),
                              attr(locfit_aggr[[1]][[1]],"prev15to49")),
                     model = "eppasm",
                     variable = c(rep("Incidence",times = (2021-1970)+1),
                                    rep("Prevalence",times = (2021-1970)+1)))

nat_dat = data.table(year = 1970:2021,
                     value= c(g1,g2),
                     model = "eppasm - with turnover",
                     variable = c(rep("Incidence",times = (2021-1970)+1),
                                  rep("Prevalence",times = (2021-1970)+1)))
epp[variable == "prev", variable := "Prevalence"]
epp[variable == "incid", variable := "Incidence"]
compare_all = rbind(nat_dat, epp[variable %in% c("Prevalence","Incidence")])
compare_all = rbind(prev_dat, compare[variable == "Prevalence",.(year, model,prev = value, type = "estimate")], fill=TRUE)

compare_all = compare_all[year <= 2020 & year >= 1994]
ggplot() + 
  geom_line(data = compare_all[type == "estimate"],aes(year,prev,col=model )) +
  geom_point(data = compare_all[type == "anc"],aes(year,prev), alpha=0.4) + 
  geom_point(data = compare_all[type == "hhs"],aes(year,prev)) + 
  geom_errorbar(data = compare_all[type == "hhs"],aes(x = year, ymin = prev - 1.96*se,
                                                      ymax = prev + 1.96*se)) +
  theme_bw() + ylab("Prevalence Rate")

ggplot() + geom_line(data = compare_all,aes(year,value,col=model),size = 1) + 
  theme_bw() + facet_wrap(~variable,scales="free_y") + ylab("Rate") + 
  scale_color_manual(values = c("maroon","turquoise"))

##With turnover
t1 = readRDS("pre-adjust.RDS")
#t1 = get_dat(fit1[[group]],"EPPASM-No Turnover")[variable == "prev"]
t2 = locfit_aggr$allmod[[group]]
prev1 = attr(t2[[1]],"hivpop")
prev1 = apply(prev1[,1:8,,],3:4,sum)[1,]

prev2 = apply(t1[,1:8,,],3:4,sum)[1,]

age15to49 = attr(t2[[1]],"allpop")[locfit$`Pop féminine restante`$fp$ss$p.age15to49.idx,,,]
popa = apply(age15to49, 4, sum)
pop1 = data.table(value = prev1/popa, variable = "prev", model = "EPPASM-Turnover",year = 1:52)
pop2 = data.table(value = prev2/popa, variable = "prev", model = "EPPASM-No Turnover",year = 1:52)

combine = rbind(pop1, pop2)
ggplot(combine, aes(year, value, col=model)) + geom_line() + theme_bw() + ylab("Prevalence Rate")

