

rm(list=ls())
devtools::install_github("mrc-ide/eppasm@new-master-csvar")
devtools::install_github("jeffeaton/epp@deepa-dev")
devtools::load_all("~/Documents/eppasm")

pjnz <- system.file("extdata/testpjnz", "senegal_2021_05_04.pjnz", package="eppasm")
entry_rates = read.csv(system.file("extdata/entry_rates.csv", package = "eppasm"))
pops = epp::read_spt(system.file("extdata/testpjnz", "senegal_2021_05_04.pjnz", package="eppasm"))
entrantprev = readRDS(system.file("extdata/entrantprev.RDS", package="eppasm"))
sen <- prepare_spec_fit(pjnz, proj.end=2019.5, popupdate = FALSE)
attr(sen$`Travailleurs du sexe`,"specfp")$popadjust <- FALSE

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

group = "Travailleurs du sexe"
nn1 = readRDS("~/Documents/tests/SEN_test8fitmod_all.RDS")
nn2 = readRDS("~/Documents/tests/SEN_test9fitmod_all.RDS")

get_mod <- function(fit){
  
  if("resample" %in% names(fit)){
    theta_ur <- apply(fit$resample,2,mean)
  } else {
    theta_ur <- fit$par
  }
  
  param <- fnCreateParam(theta_ur, fit$fp)
  fp_par <- update(fit$fp, list = param)
  
  mod <- simmod.specfp(fp_par, VERSION = "R")

  return(mod)
}

## Age distribution plots
plot_age_dist <- function(mod, lab, title, sex_both=FALSE, plhiv = FALSE){

  if(plhiv){
    ff = alply(attr(mod,population)[,,2,year1],3)
    ff = lapply(ff, function(x) prop.table(x))
  } else {
    ff = alply(attr(mod,population)[,,,year1],4)
    ff = lapply(ff, function(x) prop.table(apply(x,1:2,sum)))
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
    geom_bar(stat="identity",position = "dodge",alpha=0.8) + 
    theme_bw() + ylab(lab) + xlab("Age") + ggtitle(title)  
  
  return(gg)
  
}

fit1 = nn1[[group]]
fit2 = nn2[[group]]

population = "allpop"
plhiv = T
lab = "proportion of population"
title = "Age distribution of PLHIV - FSW"
gg1 = plot_age_dist(get_mod(fit1), lab = lab, title = title, plhiv)

title = "Age distribution of PLHIV - FSW, no ANC data"
gg2 = plot_age_dist(get_mod(fit2), lab = lab, title = title, plhiv)

##Compare to data & EPP results
group = "Clients des travailleurs du sexe"
fit1 = readRDS("~/Documents/tests/SEN_test10fitmod_all.RDS")
fit2 = readRDS("~/Documents/tests/SEN_test10fitmod_all.RDS")
pop = epp::read_spt("~/Documents/HIV/data/pjnz/senegal_2021_05_04.pjnz")
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
}

compare = rbind(get_dat(fit1[[group]],"eppasm"), 
                #get_dat(fit2[[group]],"eppasm-noANC"),
                epp,fill=TRUE,use.names=TRUE)

ggplot(compare[year <= 48], aes(year, value, col = model)) + geom_line() +
  facet_wrap(~variable, scales = "free_y") + theme_bw()
anc = data.table(attr(sen[[group]],"eppd")$ancsitedat)[,.(year,prev,type="anc",se=NA)]
hhs = data.table(attr(sen[[group]],"eppd")$hhs)[,.(year,prev,se,type="hhs")]
prev_dat = rbind(anc,hhs)

compare[,year := year + 1969]

compare_all = rbind(prev_dat, compare[variable == "prev",.(year, model,prev = value, type = "estimate")], fill=TRUE)
compare_all = compare_all[year <= 2017]
ggplot() + 
  geom_line(data = compare_all[type == "estimate"],aes(year,prev,col=model )) +
  geom_point(data = compare_all[type == "anc"],aes(year,prev), alpha=0.4) + 
  geom_point(data = compare_all[type == "hhs"],aes(year,prev)) + 
  geom_errorbar(data = compare_all[type == "hhs"],aes(x = year, ymin = prev - 1.96*se,
                                                      ymax = prev + 1.96*se)) +
  theme_bw()



