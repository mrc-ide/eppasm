setwd("~/Documents/Research/all-cause-mortality/spectrum/")

########################
####  GAM analysis  ####
########################

## Load processed DHS sibling history data
sib.mx.tips <- readRDS("~/Documents/Research/all-cause-mortality/sibling-analysis/sib-mx-tips_2017-03-09.RData")


## Define TIPS as factor

sibmx <- subset(sib.mx.tips, country %in% c("Zimbabwe", "Malawi", "Kenya", "Madagascar", "Mali", "Namibia", "Rwanda", "Tanzania", "Uganda", "Zambia", "Burkina Faso", "Cameroon", "Chad", "Cote d'Ivoire", "Ethiopia", "Guinea", "Lesotho", "Mozambique", "Niger", "Senegal", "Benin", "Congo", "DRC", "Gabon", "Liberia", "Nigeria", "Sierra Leone", "Togo"))
sibmx$ftips <- factor(sibmx$tips)
sibmx$country <- factor(sibmx$country)

library(mgcv)
library(parallel)
options(mc.cores=6)

mod.nb.k8 <- mclapply(split(sibmx, factor(sibmx$country)),
                      function(dat) gam(deaths ~ -1 + ftips + sex + s(agegr, by=sex, k=8) + s(period, by=sex, k=8) + ti(agegr, period, by=sex, k=c(8,8)), nb, dat, offset=log(pys)))

mod.nb.ka13t8 <- mclapply(split(sibmx, factor(sibmx$country)),
                          function(dat) gam(deaths ~ -1 + ftips + sex + s(agegr, by=sex, k=13) + s(period, by=sex, k=8) + ti(agegr, period, by=sex, k=c(13,8)), nb, dat, offset=log(pys)))

save(mod.nb.k8, file="gam-nb-k8-fits_2016-03-09.RData")
save(mod.nb.ka13t8, file="mod-nb-ka13t8-fits_2016-03-12.RData")

load("gam-nb-k8-fits_2016-03-09.RData")
load("mod-nb-ka13t8-fits_2016-03-12.RData")

coef.k8 <- lapply(mod.nb.k8, coef)
coef.k8 <- mapply("[", coef.k8, lapply(lapply(coef.k8, names), grep, pattern="ftips"), SIMPLIFY=FALSE)

coef.ka13t8 <- lapply(mod.nb.ka13t8, coef)
coef.ka13t8 <- mapply("[", coef.ka13t8, lapply(lapply(coef.ka13t8, names), grep, pattern="ftips"), SIMPLIFY=FALSE)

tips.mat <- do.call(rbind, coef.k8)
names(dimnames(tips.mat)) <- c("country", "tips")

Atips <- t(diag(15) - rep(c(0, 1/3, 0), c(1, 3, 11)))
dimnames(Atips) <- list(tips=0:14, tips=0:14)

tips.coef <- as.data.frame.table(Atips %*% t(tips.mat), responseName="tipscoef")

sibmx <- within(merge(sibmx, tips.coef), adjdeaths <- deaths/exp(tipscoef))

adjagemx <- within(aggregate(cbind(adjdeaths, pys) ~ period+agegr+sex+country, sibmx, sum), adjmx <- adjdeaths/pys)
adjq4515 <- within(aggregate(cbind(cummx4515=ifelse(agegr %in% 15:59, adjmx, 0), cummx3515=ifelse(agegr %in% 15:49, adjmx, 0)) ~ period+sex+country, adjagemx, sum),{q3515 <- 1-exp(-cummx3515); q4515 <- 1-exp(-cummx4515)})

devtools::use_data(adjagemx, adjq4515, pkg="~/Documents/Code/R/eppspectrum/", overwrite=TRUE)


## Calculate predicted 45q15 trend

pred_nqx<- function(mod, ..., agevar="agegr", agerange=15:59, tipscomb=NULL){

  ## @tipscomb: linear combination of TIPS for reference category, if not NULL

  newdata <- merge(data.frame(...), setNames(data.frame(agerange), agevar))
  X <- predict(mod, newdata, type="lpmatrix")
  if(!is.null(tipscomb))
    X[,grep("ftips", colnames(X))] <- rep(tipscomb, each=nrow(X))
  mx <- exp(as.numeric(X %*% coef(mod)))
  se <- sqrt(sum(diag(mx) %*% X %*% vcov(mod) %*% t(X) %*% diag(mx)))

  return(setNames(1-exp(sum(-mx) + c(0, 1, -1)*se*qnorm(0.975)),
                  c(paste0("q", diff(range(agerange))+1, min(agerange)),
                    "ci_l", "ci_u")))
}

trend_nqx <- function(mod, tipsref=1, agerange=15:59, tipscomb=rep(c(0, 1/3, 0), c(1, 3, 11))){


  newdata <-  merge(unique(mod$model[c("sex", "period")]), data.frame(ftips=as.character(tipsref)))
  newdata <- with(newdata, newdata[order(sex, period, ftips),])
  tmp <- data.frame(newdata, t(do.call(mapply, c(list(FUN=pred_nqx, MoreArgs=list(mod=mod, agerange=agerange, tipscomb=tipscomb)), newdata))))
}


sibq4515 <- lapply(mod.nb.k8, trend_nqx)
sibq4515 <- do.call(rbind, mapply("[[<-", sibq4515, "country", value=names(sibq4515), SIMPLIFY=FALSE))
sibq4515$country <- factor(sibq4515$country)

## sibq4515_ka13t8 <- mclapply(mod.nb.ka13t8, trend_nqx)
## sibq4515_ka13t8 <- do.call(rbind, mapply("[[<-", sibq4515_ka13t8, "country", value=names(sibq4515_ka13t8), SIMPLIFY=FALSE))
## sibq4515_ka13t8$country <- factor(sibq4515_ka13t8$country)

sibq3515 <- lapply(mod.nb.k8, trend_nqx, agerange=15:49)
sibq3515 <- do.call(rbind, mapply("[[<-", sibq3515, "country", value=names(sibq3515), SIMPLIFY=FALSE))
sibq3515$country <- factor(sibq3515$country)

## sibq3515ref2 <- mclapply(mod.nb.k8, trend_nqx, agerange=15:49, tipsref=2, tipscomb=NULL)
## sibq3515ref2 <- do.call(rbind, mapply("[[<-", sibq3515ref2, "country", value=names(sibq3515ref2), SIMPLIFY=FALSE))
## sibq3515ref2$country <- factor(sibq3515ref2$country)

sibq1010 <- lapply(mod.nb.k8, trend_nqx, agerange=10:19)
sibq1010 <- do.call(rbind, mapply("[[<-", sibq1010, "country", value=names(sibq1010), SIMPLIFY=FALSE))
sibq1010$country <- factor(sibq1010$country)

sibq1040 <- lapply(mod.nb.k8, trend_nqx, agerange=40:49)
sibq1040 <- do.call(rbind, mapply("[[<-", sibq1040, "country", value=names(sibq1040), SIMPLIFY=FALSE))
sibq1040$country <- factor(sibq1040$country)


## Predicted mortality rate for ages 15:60 in years 1998, 2005, 2012, for TIPS=2

predagemx <- function(mod, tipsref=2){
  newdata <- expand.grid(sex=c("male", "female"), period=unique(mod$model$period), ftips=as.character(tipsref), agegr=5:70)
  agemxpred <- setNames(predict(mod, newdata, se.fit=TRUE), c("lmx", "se.lmx"))
  data.frame(newdata, agemxpred)
}

sibagemx <- lapply(mod.nb.k8, predagemx)
sibagemx <- do.call(rbind, mapply("[[<-", sibagemx, "country", value=names(sibagemx), SIMPLIFY=FALSE))
sibagemx$country <- factor(sibagemx$country)


devtools::use_data(sibq4515, sibq3515, sibq1010, sibq1040, sibagemx, pkg="~/Documents/Code/R/eppspectrum/", overwrite=TRUE)




#######################################
####  DHS recent household deaths  ####
#######################################

library(survey)

hhd <- readRDS("~/Documents/Research/all-cause-mortality/sibling-analysis/dhs-hhdeaths_2017-03-11.RDS")
hhd$year <- as.integer(substr(hhd$survyear, 1, 4))

## Design based estiamtes
hhd_des <- svydesign(~hv002, strata=~country+survyear+year+stratum, nest=TRUE, data=hhd, weights=~weight)

hhd_mx <- svyby(~death, ~country+sex+survyear+year+age, subset(hhd_des, age %in% 15:59), svyratio, denominator=~pys, vartype="var")
names(hhd_mx)[names(hhd_mx)=="death/pys"] <- "mx"
names(hhd_mx)[names(hhd_mx)=="var.death/pys"] <- "var"

hhd_q4515 <- aggregate(cbind(mx, var) ~ country+sex+survyear+year, hhd_mx, sum)
hhd_q4515 <- cbind(hhd_q4515, 1-exp(-with(hhd_q4515, mx + qnorm(0.975)*sqrt(var) %*% t(c(q4515=0, ci_l=-1, ci_u=1)))))

hhd_q3515 <- aggregate(cbind(mx, var) ~ country+sex+survyear+year, subset(hhd_mx, age %in% 15:49), sum)
hhd_q3515 <- cbind(hhd_q3515, 1-exp(-with(hhd_q3515, mx + qnorm(0.975)*sqrt(var) %*% t(c(q3515=0, ci_l=-1, ci_u=1)))))

## normalize weights

totweight <- aggregate(cbind(totweight=weight, totpys=pys) ~ country+survyear, hhd, sum)
totweight # close to normalized already, use as is

## Smoothed age-specific mortality
hhd$age <- as.integer(as.character(hhd$age))

hhdaggr <- aggregate(weight*cbind(death, pys) ~ country+survyear+year+sex+age, hhd, sum)



hhdmod.nb.k13 <- lapply(split(hhdaggr, factor(with(hhdaggr, paste(country, year, sex, sep=".")))),
                        function(dat) gam(death ~ s(age, k=13), nb, dat, offset=log(pys)))


hhd_q4515_gam <- mapply(pred_nqx, hhdmod.nb.k13,
                        country=sub("([^/.]*).([^/.]*).([^/.]*)", "\\1", names(hhdmod.nb.k13)),
                        year=sub("([^/.]*).([^/.]*).([^/.]*)", "\\2", names(hhdmod.nb.k13)),
                        sex=sub("([^/.]*).([^/.]*).([^/.]*)", "\\3", names(hhdmod.nb.k13)),
                        agevar="age")

hhd_q4515_gam <- data.frame(country=sub("([^/.]*).([^/.]*).([^/.]*)", "\\1", colnames(hhd_q4515_gam)),
                            year=sub("([^/.]*).([^/.]*).([^/.]*)", "\\2", colnames(hhd_q4515_gam)),
                            sex=sub("([^/.]*).([^/.]*).([^/.]*)", "\\3", colnames(hhd_q4515_gam)),
                            t(hhd_q4515_gam))

merge(hhd_q4515[c("country", "sex", "survyear", "year", "q4515", "ci_l", "ci_u")], hhd_q4515_gam,
      by=c("country", "sex", "year"), suffixes=c("", "_gam"))

hhd_agemx <- mapply(function(mod, newdata, age=5:70) data.frame(newdata, setNames(predict(mod, newdata, se.fit=TRUE), c("lmx", "se.lmx"))), hhdmod.nb.k13,
                    lapply(mapply(data.frame, country=sub("([^/.]*).([^/.]*).([^/.]*)", "\\1", names(hhdmod.nb.k13)),
                                  year=as.integer(sub("([^/.]*).([^/.]*).([^/.]*)", "\\2", names(hhdmod.nb.k13))),
                                  sex=sub("([^/.]*).([^/.]*).([^/.]*)", "\\3", names(hhdmod.nb.k13)), SIMPLIFY=FALSE),
                           merge, data.frame(age=5:70)),
                    SIMPLIFY=FALSE)
hhd_agemx <- plyr::rbind.fill(hhd_agemx)


######################################
####  Load UNAIDS 2016 estimates  ####
######################################

devtools::load_all("~/Documents/Code/R/eppspectrum/")

specdir <- "~/Documents/Data/Spectrum files/2016 final/SSA"

specres <- list(Benin = aggr_specres(lapply(list.files(specdir, "Benin", full.names=TRUE), read_hivproj_output)),
                "Burkina Faso" = read_hivproj_output(list.files(specdir, "Burkina", full.names=TRUE)),
                Cameroon = read_hivproj_output(list.files(specdir, "Cameroon", full.names=TRUE)),
                Chad = read_hivproj_output(list.files(specdir, "Chad", full.names=TRUE)),
                Congo = read_hivproj_output(list.files(specdir, "Congo", full.names=TRUE)),
                "Cote d'Ivoire" = aggr_specres(lapply(list.files(specdir, "Cote", full.names=TRUE), read_hivproj_output)),
                DRC = read_hivproj_output(list.files(specdir, "DRC", full.names=TRUE)),
                Ethiopia      = read_hivproj_output(list.files(specdir, "Ethiopia", full.names=TRUE)),
                Gabon         = read_hivproj_output(list.files(specdir, "Gabon", full.names=TRUE)),
                Guinea        = read_hivproj_output(list.files(specdir, "Guinea", full.names=TRUE)),
                Kenya         = aggr_specres(lapply(list.files(specdir, "Kenya", full.names=TRUE), read_hivproj_output)),
                Lesotho       = read_hivproj_output(list.files(specdir, "Lesotho", full.names=TRUE)),
                Liberia       = read_hivproj_output(list.files(specdir, "Liberia", full.names=TRUE)),
                ## Madagascar    = read_hivproj_output(list.files(specdir, "Madagascar", full.names=TRUE)),
                Malawi        = read_hivproj_output(list.files(specdir, "Malawi", full.names=TRUE)[1]),
                Mali          = read_hivproj_output(list.files(specdir, "Mali", full.names=TRUE)),
                ## Mozambique    = aggr_specres(lapply(list.files(specdir, "Moz", full.names=TRUE), read_hivproj_output)),
                Namibia       = read_hivproj_output(list.files(specdir, "Namibia", full.names=TRUE)),
                Niger         = read_hivproj_output(list.files(specdir, "^Niger_", full.names=TRUE)),
                Nigeria       = aggr_specres(lapply(list.files(specdir, "Nigeria", full.names=TRUE), read_hivproj_output)),
                Rwanda        = read_hivproj_output(list.files(specdir, "Rwanda", full.names=TRUE)),
                Senegal       = read_hivproj_output(list.files(specdir, "Senegal", full.names=TRUE)),
                "Sierra Leone"= read_hivproj_output(list.files(specdir, "Sierra", full.names=TRUE)),
                Tanzania      = read_hivproj_output(list.files(specdir, "TZ", full.names=TRUE)),
                Togo          = read_hivproj_output(list.files(specdir, "Togo", full.names=TRUE)),
                Uganda        = read_hivproj_output(list.files(specdir, "Uganda", full.names=TRUE)[1]),
                Zambia        = read_hivproj_output(list.files(specdir, "Zambia", full.names=TRUE)),
                Zimbabwe      = aggr_specres(lapply(list.files(specdir, "Zimbabwe", full.names=TRUE), read_hivproj_output)))


spec16q4515 <- lapply(specres, calc_nqx)
spec16natq4515 <- lapply(specres, calc_nqx, nonhiv=TRUE)

spec16q3515 <- lapply(specres, calc_nqx, n=35)
spec16natq3515 <- lapply(specres, calc_nqx, n=35, nonhiv=TRUE)

spec16agemx <- lapply(specres, agemx)

############################
####  Plot 45q15 trend  ####
############################

library(RColorBrewer)

cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[,1], rev(y[,2])), border=NA, ...)

transp <- function(col, alpha=0.5)
  return(apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha)))



## X11(h=5.5, w=12, pointsize=18, type="cairo")
## quartz(h=5.5, w=12, pointsize=18)


## Plot 45q15

plotdata <- TRUE
plotaggr <- FALSE
plotagemod <- FALSE
plotspecres <- TRUE

if(plotspecres){ 
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q4515_sibmx-unaids2016_2017-03-13.pdf"
} else { 
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q4515_sibmx_2017-03-13.pdf"
}

pdf(filename, h=5.5, w=12, pointsize=18)
##
for(icountry in levels(sibq4515$country)){
  par(mfrow=c(1, 2), las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.0, 2, 0.5), cex=0.8)
  ##
  for(isex in rev(levels(sibq4515$sex))){
    dat <- subset(sibq4515, sex==isex & country == icountry)
    cols <- brewer.pal(3, "Set1")
    plot(dat$period, dat$q4515, type="n",
         ylim=c(0, 0.9), xlim=c(1985, 2015), ylab=expression(""[45]~q[15]), xlab="", main=paste(icountry, c("male"="Men", "female"="Women")[isex]))
    if(plotdata)
      with(subset(adj45q15, sex == isex & country == icountry), points(period, q4515, pch=23, col=transp("grey", 0.4), bg=transp("grey", 0.4)))
    cred.region(dat$period, dat[c("ci_l", "ci_u")], col=transp(cols[2], 0.3))
    ## if(plotaggr){
    ##   est <- estci(q4515$aggr_nopreg[[icountry]][isex,as.character(1983:2015),])
    ##   cred.region(1983:2015, est[,c("lower", "upper")], col=transp("darkolivegreen", 0.3))
    ##   lines(1983:2015, est[,"mean"], col="darkolivegreen", lwd=2)
    ## }
    ## if(plotagemod){
    ##   est <- estci(q4515$ageprev[[icountry]][isex,as.character(1983:2015),])
    ##   cred.region(1983:2015, est[,c("lower", "upper")], col=transp("darkred", 0.3))
    ##   lines(1983:2015, est[,"mean"], col="darkred", lwd=2)
    ## }
    if(plotspecres){
      lines(1983:2015, spec16q4515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=2, lwd=3)
      lines(1983:2015, spec16natq4515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=3, lwd=3, col="grey30")
    }
    ## if(plothhd){
    ##   with(subset(hhd_q4515, country == icountry & sex == isex), segments(year, y0=ci_l, y1=ci_u, lwd=2))
    ##   with(subset(hhd_q4515, country == icountry & sex == isex), points(year, q4515, pch=19))
    ## }
    lines(dat$period, dat$q4515, lwd=2, lty=1, col=cols[2])
    ## if(plotaggr & !plotagemod)
    ##   legend("topleft", c("Sibling history", "Default IRRs"), pch=15, pt.cex=2.2, col=transp(c(cols[2], "darkolivegreen")), bty="n")
    ## if(plotagemod)
    ##   legend("topleft", c("Sibling history", "Default IRRs", "Estimated IRRs", "UNAIDS 2016"), pch=c(15, 15, 15, NA), pt.cex=2.2, col=c(transp(c(cols[2], "darkolivegreen", "darkred")), 1), lty=c(NA, NA, NA, 2), lwd=c(NA, NA, NA, 2), bty="n")
    ## legend("topleft", c("Sibling history", "Estimated IRRs", "UNAIDS 2016", "UNAIDS '16 non-HIV"), pch=c(15, 15, NA, NA), pt.cex=2.2, col=c(transp(c(cols[2], "darkred")), 1, "grey30"), lty=c(NA, NA, 2, 3), lwd=c(NA, NA, 2, 2), bty="n", cex=0.85)
  }
}
##
dev.off()



## Plot 35q15

plotdata <- TRUE
plotaggr <- FALSE
plotagemod <- FALSE
plotspecres <- FALSE

if(plotspecres){ 
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q3515_sibmx-unaids2016_2017-03-13.pdf"
} else { 
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q3515_sibmx_2017-03-13.pdf"
}

pdf(filename, h=5.5, w=12, pointsize=18)
##
for(icountry in levels(sibq3515$country)){
  par(mfrow=c(1, 2), las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.0, 2, 0.5), cex=0.8)
  ##
  for(isex in rev(levels(sibq3515$sex))){
    dat <- subset(sibq3515, sex==isex & country == icountry)
    cols <- brewer.pal(3, "Set1")
    plot(dat$period, dat$q3515, type="n",
         ylim=c(0, 0.8), xlim=c(1985, 2015), ylab=expression(""[35]~q[15]), xlab="", main=paste(icountry, c("male"="Men", "female"="Women")[isex]))
    if(plotdata)
      with(subset(adj45q15, sex == isex & country == icountry), points(period, q3515, pch=23, col=transp("grey", 0.4), bg=transp("grey", 0.4)))
    cred.region(dat$period, dat[c("ci_l", "ci_u")], col=transp(cols[2], 0.3))
    ## if(plotaggr){
    ##   est <- estci(q3515$aggr_nopreg[[icountry]][isex,as.character(1983:2015),])
    ##   cred.region(1983:2015, est[,c("lower", "upper")], col=transp("darkolivegreen", 0.3))
    ##   lines(1983:2015, est[,"mean"], col="darkolivegreen", lwd=2)
    ## }
    if(plotagemod){
      est <- estci(q3515$ageprev[[icountry]][isex,as.character(1983:2015),])
      cred.region(1983:2015, est[,c("lower", "upper")], col=transp("darkred", 0.3))
      lines(1983:2015, est[,"mean"], col="darkred", lwd=2)
    }
    if(plotspecres){
      lines(1983:2015, spec16q3515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=2, lwd=3)
      lines(1983:2015, spec16natq3515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=3, lwd=3, col="grey30")
    }
    ## if(plothhd){
    ##   with(subset(hhd_q3515, country == icountry & sex == isex), segments(year, y0=ci_l, y1=ci_u, lwd=2))
    ##   with(subset(hhd_q3515, country == icountry & sex == isex), points(year, q3515, pch=19))
    ## }
    lines(dat$period, dat$q3515, lwd=2, lty=1, col=cols[2])
    ## if(plotaggr & !plotagemod)
    ##   legend("topleft", c("Sibling history", "Default IRRs"), pch=15, pt.cex=2.2, col=transp(c(cols[2], "darkolivegreen")), bty="n")
    ## if(plotagemod)
    ##   ## legend("topleft", c("Sibling history", "Default IRRs", "Estimated IRRs"), pch=15, pt.cex=2.2, col=transp(c(cols[2], "darkolivegreen", "darkred")), bty="n")
    ## legend("topleft", c("Sibling history", "Estimated IRRs", "UNAIDS 2016", "UNAIDS '16 non-HIV"), pch=c(15, 15, NA, NA), pt.cex=2.2, col=c(transp(c(cols[2], "darkred")), 1, "grey30"), lty=c(NA, NA, 2, 3), lwd=c(NA, NA, 2, 2), bty="n", cex=0.85)
  }
}
##
dev.off()


########################################
####  Household deaths comparison   ####
########################################

## X11(h=5, w=12, pointsize=14, type="cairo")

## 35q15
if(plotspecres){
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q3515_sibmx-hhdeath-specres-compare_2017-03-13.pdf"
 } else {
   filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q3515_sibmx-hhdeath-compare_2017-03-13.pdf"
 }
pdf(filename, h=5, w=12, pointsize=14)
##
par(mfrow=c(2, 6), las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, oma=c(0, 3, 0, 0), mar=c(2.0, 0.5, 1.5, 0.5), cex=0.7)
layout(rbind(1:6, 7:12, 13), heights=c(1, 1, 0.2))
for(icountry in c("Cote d'Ivoire", "Nigeria", "Tanzania", "Uganda", "Malawi", "Zimbabwe")){
  ##
  for(isex in rev(levels(sibq3515$sex))){
    dat <- subset(sibq3515, sex==isex & country == icountry)
    cols <- brewer.pal(3, "Set1")
    plot(dat$period, dat$q3515, type="n",
         ylim=c(0, 0.8), xlim=c(1985, 2015), ylab="", xlab="", main=paste(icountry, c("male"="Men", "female"="Women")[isex]),
         yaxt="n")
    if(icountry %in% c("Cote d'Ivoire", "Uganda") && isex == "female"){
      axis(2, xpd=NA)
      mtext(expression(""[35]~q[15]), 2, 2, xpd=NA, las=0, cex=0.7 )
    } else
      axis(2, labels=FALSE)
    if(plotdata)
      with(subset(adj45q15, sex == isex & country == icountry), points(period, q3515, pch=23, col=transp("grey", 0.3), bg=transp("grey", 0.3)))
    if(plotspecres){
      lines(1983:2015, spec16q3515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=2, lwd=2)
      lines(1983:2015, spec16natq3515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=3, lwd=2, col="grey30")
    }
    cred.region(dat$period, dat[c("ci_l", "ci_u")], col=transp(cols[2], 0.3))
    ##
    with(subset(hhd_q3515, country == icountry & sex == isex), segments(year, y0=ci_l, y1=ci_u, lwd=2))
    with(subset(hhd_q3515, country == icountry & sex == isex), points(year, q3515, pch=16))
    ##
    lines(dat$period, dat$q3515, lwd=2, lty=1, col=cols[2])
  }
}
####
par(mar=c(0, 0, 0, 0))
plot(0, 0, type="n", bty="n", axes=FALSE)
legend("center", c("Adjusted sibling history estamates   ", "Recent household deaths"), pch=c(15, 16), pt.cex=c(3, 1.0), lty=c(NA, 1), lwd=2, col=c(transp(cols[2], 0.3), 1), horiz=TRUE, bty="n")
####
dev.off()


## 45q15
plotdata <- TRUE

if(plotspecres){
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q4515_sibmx-hhdeath-specres-compare_2017-03-13.pdf"
 } else {
   filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/q4515_sibmx-hhdeath-compare_2017-03-13.pdf"
 }
pdf(filename, h=5, w=12, pointsize=14)
par(mfrow=c(2, 6), las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, oma=c(0, 3, 0, 0), mar=c(2.0, 0.5, 1.5, 0.5), cex=0.7)
layout(rbind(1:6, 7:12, 13), heights=c(1, 1, 0.2))
for(icountry in c("Cote d'Ivoire", "Nigeria", "Tanzania", "Uganda", "Malawi", "Zimbabwe")){
  ##
  for(isex in rev(levels(sibq4515$sex))){
    dat <- subset(sibq4515, sex==isex & country == icountry)
    cols <- brewer.pal(3, "Set1")
    plot(dat$period, dat$q4515, type="n",
         ylim=c(0, 0.9), xlim=c(1985, 2015), ylab="", xlab="", main=paste(icountry, c("male"="Men", "female"="Women")[isex]),
         yaxt="n")
    if(icountry %in% c("Cote d'Ivoire", "Uganda") && isex == "female"){
      axis(2, xpd=NA)
      mtext(expression(""[45]~q[15]), 2, 2, xpd=NA, las=0, cex=0.7 )
    } else
      axis(2, labels=FALSE)
    if(plotdata)
      with(subset(adj45q15, sex == isex & country == icountry), points(period, q4515, pch=23, col=transp("grey", 0.3), bg=transp("grey", 0.3)))
    cred.region(dat$period, dat[c("ci_l", "ci_u")], col=transp(cols[2], 0.3))
    if(plotspecres){
      lines(1983:2015, spec16q4515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=2, lwd=2)
      lines(1983:2015, spec16natq4515[[icountry]][c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=3, lwd=2, col="grey30")
    }
    ##
    with(subset(hhd_q4515, country == icountry & sex == isex), segments(year, y0=ci_l, y1=ci_u, lwd=2))
    with(subset(hhd_q4515, country == icountry & sex == isex), points(year, q4515, pch=16))
    ##
    lines(dat$period, dat$q4515, lwd=2, lty=1, col=cols[2])
  }
}
####
par(mar=c(0, 0, 0, 0))
plot(0, 0, type="n", bty="n", axes=FALSE)
legend("center", c("Adjusted sibling history estamates   ", "Recent household deaths"), pch=c(15, 16), pt.cex=c(3, 1.0), lty=c(NA, 1), lwd=2, col=c(transp(cols[2], 0.3), 1), horiz=TRUE, bty="n")
####
dev.off()


## Compare agemx

agemx_compare <- merge(hhd_agemx, sibagemx, by.x=c("country", "year", "sex", "age"), by.y=c("country", "period", "sex", "agegr"), suffixes=c("_hhd", "_sib"))
groups <- unique(agemx_compare[c("country", "year", "sex")])
groups <- subset(groups, country %in% c("Malawi", "Tanzania", "Uganda", "Zimbabwe"))


## pdf("figures/sib-hhd-agemx-compare_2017-03-12.pdf", h = 6.2, w = 12, pointsize = 18)
pdf("~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/sib-hhd-unaids16_agemx-compare_2017-03-13.pdf", h = 6.2, w = 12, pointsize = 18)
##
par(las=1, mgp=c(1.5, 0.5, 0), tcl=-0.25, oma=c(0, 3.0, 0, 0.5), mar=c(0.5, 0.5, 0.25, 0.25), cex=0.8)
layout(rbind(1:4, 5:8, 0, 9), heights=c(1,1,0.1, 0.2))
for(ii in 1:nrow(groups)){
  icountry <- as.character(groups$country[ii])
  iyear <- groups$year[ii]
  isex <- groups$sex[ii]
  ##
  dat <- subset(agemx_compare, sex == isex & year == iyear & age %in% 15:59 & country == icountry)
  with(dat, plot(age, lmx_hhd, type="n", ylim=c(-7, -2.5), ylab="", xlab="", main="", xaxt="n", yaxt="n"))
  mtext(paste(icountry, c(male="Men", female="Women")[isex], iyear), 3, -1.1, adj=0.05, font=2, cex=0.8)
  yaxl <- c(1, 2, 5, 10, 25, 50)
  yax <- log(yaxl/1000)
  axis(1, labels=FALSE); axis(2, yax, labels=FALSE)
  if((ii-1) %% 8 >= 4 || icountry == "Uganda"){axis(1, tick=FALSE); mtext("Age", 1, 1.5, cex=0.8)}
  if(ii %% 4 == 1){axis(2, yax, yaxl, tick=FALSE)}
  with(dat, cred.region(age, c(lmx_hhd) + (qnorm(0.975)*se.lmx_hhd %*% t(c(-1,1))), col=transp("darkred", 0.3)))
  with(dat, cred.region(age, c(lmx_sib) + (qnorm(0.975)*se.lmx_sib %*% t(c(-1,1))), col=transp(cols[2], 0.3)))
  if(plotspecres){
    if(!is.null(spec16agemx[[icountry]]))
      lines(15:59, log(spec16agemx[[icountry]][as.character(15:59), c("male"="Male", "female"="Female")[isex], as.character(iyear)]), lty=2, lwd=2)
    ## lines(15:49, log(spec16natagemx[[icountry]][as.character(15:49), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=3, lwd=2, col="grey30")
  }
  with(dat, matlines(age, cbind(lmx_hhd, lmx_sib), col=c("darkred", cols[2]), lwd=2, lty=1))
  if(ii == 14){
    plot(0, 0, type="n", bty="n", axes=FALSE, xlab="", ylab="")
    plot(0, 0, type="n", bty="n", axes=FALSE, xlab="", ylab="")
    ii <- 16
  }
  if(ii %% 8 == 0){
    par(mar=rep(0, 4))
    plot(0, 0, type="n", bty="n", axes=FALSE)
    legend("bottom", c("Adjusted sibling history", "Recent household deaths"), pch=c(15, 15), pt.cex=2.5, col=transp(c(cols[2], "darkred")), horiz=TRUE, bty="n")
    mtext("Deaths per 1000 (log scale)", 2, 1.5, cex=0.8, las=0, outer=TRUE)
    par(mar=c(0.5, 0.5, 0.25, 0.25))
  }
}
##
dev.off()



#######################
####  Plot age mx  ####
#######################

periods <- round(mapply(seq,
                        from=pmax(with(sibagemx, tapply(period, country, min)+2), 1988),
                        to=pmin(with(sibagemx, tapply(period, country, max)), 2015),
                        length.out=4))
periods[4,"Rwanda"] <- 2014

plotaggr <- FALSE
plotagemod <- FALSE


if(plotspecres){
  filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/agemx_sibmx-unaids2016_2017-03-13.pdf"
 } else {
   filename <- "~/Documents/Meetings/2017/WHO Global Health Statistics reference group/analysis/agemx_sibmx_2017-03-13.pdf"
 }
pdf(filename, h=6.5, w=12, pointsize=18)
##
for(icountry in levels(sibq4515$country)){
  par(mfrow=c(2, 4), las=1, mgp=c(1.5, 0.5, 0), tcl=-0.25, oma=c(0, 3.0, 2.2, 0.5), mar=c(0.5, 0.5, 0.25, 0.25), cex=0.8)
  layout(rbind(1:4, 5:8, 0, 9), h=c(1,1,0.1, 0.2))
  ##
  for(isex in rev(levels(sibagemx$sex))){
    iper <- periods[,icountry]
    for(per in iper){
      dat <- subset(sibagemx, sex == isex & period == per & agegr %in% 15:49 & country == icountry)
      with(dat, plot(agegr, lmx, type="n", ylim=c(-7, -2.5), ylab="", xlab="", main="", xaxt="n", yaxt="n"))
      mtext(paste(c(male="Men", female="Women")[isex], per), 3, -1.1, adj=0.05, font=2, cex=0.8)
      yaxl <- c(1, 2, 5, 10, 25, 50)
      yax <- log(yaxl/1000)
      axis(1, labels=FALSE); axis(2, yax, labels=FALSE)
      if(isex == "male"){axis(1, tick=FALSE); mtext("Age", 1, 1.5, cex=0.8)}
      if(per == iper[1]){axis(2, yax, yaxl, tick=FALSE)}
      with(dat, cred.region(agegr, c(lmx) + (qnorm(0.975)*se.lmx %*% t(c(-1,1))), col=transp(cols[2], 0.3)))
      ## if(plotaggr){
      ##   logmx <- estci(log(agemx$aggr_nopreg[[icountry]][, isex, as.character(per),]))
      ##   est <- apply(logmx, 2, function(y) predict(smooth.spline(15:80, y, spar=0.5), 15:49)$y)
      ##   cred.region(15:49, est[,c("lower", "upper")], col=transp("darkolivegreen", 0.3))
      ##   lines(15:49, est[,"mean"], col="darkolivegreen", lwd=2)
      ## }
      ## if(plotagemod){
      ##   logmx <- estci(log(agemx$ageprev[[icountry]][, isex, as.character(per),]))
      ##   est <- apply(logmx, 2, function(y) predict(smooth.spline(15:80, y, spar=0.5), 15:49)$y)
      ##   cred.region(15:49, est[,c("lower", "upper")], col=transp("darkred", 0.3))
      ##   lines(15:49, est[,"mean"], col="darkred", lwd=2)
      ## }
      if(plotspecres){
        if(!is.null(spec16agemx[[icountry]]))
          lines(15:49, log(spec16agemx[[icountry]][as.character(15:49), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=2, lwd=2)
        ## lines(15:49, log(spec16natagemx[[icountry]][as.character(15:49), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=3, lwd=2, col="grey30")
      }
      with(dat, lines(agegr, lmx, col=cols[2], lwd=2, lty=1))
    }
  }
  par(mar=rep(0, 4))
  plot(0, 0, type="n", bty="n", axes=FALSE)
  ## if(!plotaggr & !plotagemod)
  ##   legend("bottom", c("Sibling history"), pch=15, pt.cex=2.5, col=transp(c(cols[2])), horiz=TRUE, inset=0, bty="n")
  ## if(plotaggr & !plotagemod)
  ##   legend("bottom", c("Sibling history", "Default IRRs"), pch=15, pt.cex=2.5, col=transp(c(cols[2], "darkolivegreen")), horiz=TRUE, inset=0, bty="n")
  ## if(plotagemod)
  ##   ## legend("bottom", c("Sibling history", "Default IRRs", "Estimated IRRs"), pch=15, pt.cex=2.5, col=transp(c(cols[2], "darkolivegreen", "darkred")), horiz=TRUE, inset=0, bty="n")
  ## legend("bottom", c("Sibling history", "Estimated IRRs", "UNAIDS 2016", "UNAIDS '16 non-HIV"), pch=c(15, 15, NA, NA), pt.cex=2.5, col=c(transp(c(cols[2], "darkred")), 1, "grey30"), lty=c(NA, NA, 2, 3), lwd=c(NA, NA, 2, 2), bty="n", horiz=TRUE)
  mtext("Deaths per 1000 (log scale)", 2, 1.5, cex=0.8, las=0, outer=TRUE)
  mtext(icountry, 3, line=0.5, outer=TRUE, cex=1.2, font=2)
}
dev.off()



filename <- "figures/agemx_5to70_sibmx_2017-03-09.pdf"
pdf(filename, h=6.5, w=12, pointsize=18)
##
for(icountry in levels(sibq4515$country)){
  par(mfrow=c(2, 4), las=1, mgp=c(1.5, 0.5, 0), tcl=-0.25, oma=c(0, 3.0, 2.2, 0.5), mar=c(0.5, 0.5, 0.25, 0.25), cex=0.8)
  layout(rbind(1:4, 5:8, 0, 9), h=c(1,1,0.1, 0.2))
  ##
  for(isex in rev(levels(sibagemx$sex))){
    iper <- periods[,icountry]
    for(per in iper){
      dat <- subset(sibagemx, sex == isex & period == per & agegr %in% 5:70 & country == icountry)
      with(dat, plot(agegr, lmx, type="n", ylim=c(-7, -2.5), ylab="", xlab="", main="", xaxt="n", yaxt="n"))
      mtext(paste(c(male="Men", female="Women")[isex], per), 3, -1.1, adj=0.05, font=2, cex=0.8)
      yaxl <- c(1, 2, 5, 10, 25, 50)
      yax <- log(yaxl/1000)
      axis(1, labels=FALSE); axis(2, yax, labels=FALSE)
      if(isex == "male"){axis(1, tick=FALSE); mtext("Age", 1, 1.5, cex=0.8)}
      if(per == iper[1]){axis(2, yax, yaxl, tick=FALSE)}
      with(dat, cred.region(agegr, c(lmx) + (qnorm(0.975)*se.lmx %*% t(c(-1,1))), col=transp(cols[2], 0.3)))
      ## if(plotaggr){
      ##   logmx <- estci(log(agemx$aggr_nopreg[[icountry]][, isex, as.character(per),]))
      ##   est <- apply(logmx, 2, function(y) predict(smooth.spline(15:80, y, spar=0.5), 15:49)$y)
      ##   cred.region(15:49, est[,c("lower", "upper")], col=transp("darkolivegreen", 0.3))
      ##   lines(15:49, est[,"mean"], col="darkolivegreen", lwd=2)
      ## }
      if(plotagemod){
        logmx <- estci(log(agemx$ageprev[[icountry]][, isex, as.character(per),]))
        est <- apply(logmx, 2, function(y) predict(smooth.spline(15:80, y, spar=0.5), 15:49)$y)
        cred.region(15:49, est[,c("lower", "upper")], col=transp("darkred", 0.3))
        lines(15:49, est[,"mean"], col="darkred", lwd=2)
      }
      if(plotspecres){
        lines(15:49, log(spec16agemx[[icountry]][as.character(15:49), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=2, lwd=2)
        ## lines(15:49, log(spec16natagemx[[icountry]][as.character(15:49), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=3, lwd=2, col="grey30")
      }
      with(dat, lines(agegr, lmx, col=cols[2], lwd=2, lty=1))
    }
  }
  par(mar=rep(0, 4))
  plot(0, 0, type="n", bty="n", axes=FALSE)
  ## if(!plotaggr & !plotagemod)
  ##   legend("bottom", c("Sibling history"), pch=15, pt.cex=2.5, col=transp(c(cols[2])), horiz=TRUE, inset=0, bty="n")
  ## if(plotaggr & !plotagemod)
  ##   legend("bottom", c("Sibling history", "Default IRRs"), pch=15, pt.cex=2.5, col=transp(c(cols[2], "darkolivegreen")), horiz=TRUE, inset=0, bty="n")
  ## if(plotagemod)
  ##   ## legend("bottom", c("Sibling history", "Default IRRs", "Estimated IRRs"), pch=15, pt.cex=2.5, col=transp(c(cols[2], "darkolivegreen", "darkred")), horiz=TRUE, inset=0, bty="n")
  ## legend("bottom", c("Sibling history", "Estimated IRRs", "UNAIDS 2016", "UNAIDS '16 non-HIV"), pch=c(15, 15, NA, NA), pt.cex=2.5, col=c(transp(c(cols[2], "darkred")), 1, "grey30"), lty=c(NA, NA, 2, 3), lwd=c(NA, NA, 2, 2), bty="n", horiz=TRUE)
  mtext("Deaths per 1000 (log scale)", 2, 1.5, cex=0.8, las=0, outer=TRUE)
  mtext(icountry, 3, line=0.5, outer=TRUE, cex=1.2, font=2)
}
dev.off()





quartz(h=2.5, w=4, pointsize=16)

par(las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, mar=c(2, 3.0, 0.5, 0.5), cex=0.8)
matplot(1980:2020, cbind(spec16natq4515[["Zimbabwe"]]["Female", as.character(1980:2020)],
                         spec16q4515[["Zimbabwe"]]["Female", as.character(1980:2020)]), type="l", lty=3:2, lwd=2.5, col=1, ylim=c(0, 1.0), ylab=expression(""[45]~q[15]), xab="")
