plot_agefit <- function(icountry, eppmod, fitaggr, fitincrr){
  par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
  ##
  ## prevalence trend
  plot_prev(fitaggr, fitincrr, col=c("grey30", "darkred"))
  mtext("HIV prevalence, age 15-49y", 3, 0.5, font=2, cex=1.2)
  ##
  ## incidence trend
  plot_incid(fitaggr, fitincrr, col=c("grey30", "darkred"))
  mtext("HIV incidence rate, age 15-49y", 3, 0.5, font=2, cex=1.2)
  ##
  ## r-trend
  plot_rvec(fitaggr, fitincrr, col=c("grey30", "darkred"))
  mtext("r(t)", 3, 0.5, font=2, cex=1.2)
  ##
  ## vinfl
  ##
  ## sex incrr
  logincrrsex <- log(sapply(fitincrr$param, "[[", "incrr_sex")[1,])
  dens <- density(logincrrsex)
  densCI <- which(dens$x >= quantile(logincrrsex, 0.025) &
                  dens$x <= quantile(logincrrsex, 0.975))
  plot(dens, xlim=c(-0.1, 0.75), col="darkred",
       main="F:M incidence rate ratio, posterior density", xlab="log F:M IRR")
  polygon(dens$x[c(min(densCI), densCI, max(densCI))], c(0, dens$y[densCI], 0),
          col=transp("darkred"), border=NA)
  segments(x0=mean(log(sapply(fitincrr$param, "[[", "incrr_sex")[1,])), y0=0, y1=6, col="darkred", lwd=2)
  segments(x0=log(tail(fitaggr$fp$incrr_sex, 1)), y0=0, y1=6, col="grey30", lwd=2)
  lines(seq(-0.1, 0.75, 0.01), dnorm(seq(-0.1, 0.75, 0.01), eppspectrum:::sexincrr.pr.mean, eppspectrum:::sexincrr.pr.sd), col="darkblue", lty=2)
  segments(x0=eppspectrum:::sexincrr.pr.mean, y0=0, y1=2, col="darkblue", lwd=2, lty=2)
  ##
  ## age incrr
  xx <- c(1:2, 4:7)
  logincrrage <- estci(sapply(fitincrr$param, "[[", "logincrr_age"))
  ## plot men
  plot(1:7+0.5, 1:7, type="n",  xlim=c(1, 8), ylim=c(-2.5, 2),
       xlab="Age group", ylab="log IRR", main = "Male incidence rate ratios (log)", xaxt="n")
  axis(1, 1:7+0.5, paste0(3:9*5, "-", 3:9*5+4))
  abline(h=0, col="grey")
  points(3.5, 0, pch=4, lwd=2.5, col=1, cex=1.2)
  rect(xx+0.1,
       eppspectrum:::ageincrr.pr.mean[1:6]-qnorm(0.975)*eppspectrum:::ageincrr.pr.sd,
       xx+0.9,
       eppspectrum:::ageincrr.pr.mean[1:6]+qnorm(0.975)*eppspectrum:::ageincrr.pr.sd,
       col=transp("darkblue", 0.1), border=transp("darkblue", 0.5), lty=2)
  ##
  rect(xx+0.1,  logincrrage[xx,3], xx+0.9, logincrrage[xx,4], col=transp("darkred"), border=NA)
  defaultincrr <- log(fitaggr$fp$incrr_age[(xx-1)*5+1,1,fitaggr$fp$ss$PROJ_YEARS])
  segments(xx+0.1, defaultincrr, xx+0.9, col="grey30", lwd=2)
  segments(xx+0.1, eppspectrum:::ageincrr.pr.mean[1:6], xx+0.9, col=transp("darkblue", 0.5), lwd=1.5, lty=2)
  segments(xx+0.1, logincrrage[xx,1], xx+0.9, col="darkred", lwd=2)
  ## plot women
  plot(1:7+0.5, 1:7, type="n",  xlim=c(1, 8), ylim=c(-2.5, 2),
       xlab="Age group", ylab="log IRR", main = "Female incidence rate ratios (log)", xaxt="n")
  axis(1, 1:7+0.5, paste0(3:9*5, "-", 3:9*5+4))
  abline(h=0, col="grey")
  points(3.5, 0, pch=4, lwd=2.5, col=1, cex=1.2)
  rect(xx+0.1,
       eppspectrum:::ageincrr.pr.mean[7:12]-qnorm(0.975)*eppspectrum:::ageincrr.pr.sd,
       xx+0.9,
       eppspectrum:::ageincrr.pr.mean[7:12]+qnorm(0.975)*eppspectrum:::ageincrr.pr.sd,
       col=transp("darkblue", 0.1), border=transp("darkblue", 0.5), lty=2)
  ##
  rect(xx+0.1,  logincrrage[xx+7,3], xx+0.9, logincrrage[xx+7,4], col=transp("darkred"), border=NA)
  defaultincrr <- log(fitaggr$fp$incrr_age[(xx-1)*5+1,2,fitaggr$fp$ss$PROJ_YEARS])
  segments(xx+0.1, defaultincrr, xx+0.9, col="grey30", lwd=2)
  segments(xx+0.1, eppspectrum:::ageincrr.pr.mean[7:12], xx+0.9, col=transp("darkblue", 0.5), lwd=1.5, lty=2)
  segments(xx+0.1, logincrrage[xx+7,1], xx+0.9, col="darkred", lwd=2)
  ##
  mtext(paste0(icountry, ", ", eppmod, " model; posterior distribution"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  ##
  ##  Age specific prevalence compared to survey
  ##
  plot_compare_ageprev(fitaggr, fitincrr)
  mtext(paste0(icountry, ", ", eppmod, " model; Age-specific prevalence"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  ##
  ##
  ##  Age prevalence trend by age group
  ##
  par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
  ##
  for(iagegr in 1:3){
    for(isex in 1:2){
      strsex <- c("male", "female")[isex]
      stragegr <- c("15-24", "25-34", "35-49")[iagegr]
      aggr.yprev <- estci(fitaggr$agegr3prev[iagegr,isex,,])
      fit.yprev <- estci(fitincrr$agegr3prev[iagegr,isex,,])
      survdat <- subset(prev.agegr3sex.nat, country==icountry & sex == strsex & agegr3==stragegr)
      ##
      plot(1999:2016, fit.yprev[,1], type="n", ylim=c(0, 0.05*ceiling(max(survdat$ci_u, aggr.yprev[,1], fit.yprev[,1])/0.05)),
           ylab="", xlab="", main=paste(strsex, stragegr))
      cred.region(1999:2016, t(aggr.yprev[,3:4]), col=transp("grey30"))
      cred.region(1999:2016, t(fit.yprev[,3:4]), col=transp("darkred"))
      lines(1999:2016, aggr.yprev[,1], col="grey30", lwd=2)
      lines(1999:2016, fit.yprev[,1], col="darkred", lwd=2)
      ##
      points(survdat$year, survdat$prev, pch=19)
      segments(survdat$year, y0=survdat$ci_l, y1=survdat$ci_u)
    }
  }
  mtext(paste0(icountry, ", ", eppmod, " model; Prevalence trend by age group"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  ##
  ##
  return(invisible())
}
