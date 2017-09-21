plot_agefit <- function(icountry, eppmod, fitaggr, fitincrr, fit3=NULL, specres=NULL, cols=c("grey30", "darkred", "darkgreen"), pages=1:3,
                        agegr3dat = subset(prev.agegr3sex.nat, country==icountry)){
  if(1 %in% pages){
    par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
    ##
    ## prevalence trend
    if(!is.null(fit3))
      plot_prev(fitaggr, fitincrr, fit3, col=cols)
    else
      plot_prev(fitaggr, fitincrr, col=cols)
    if(!is.null(specres))
      lines(as.numeric(names(prev(specres))), prev(specres), lty=2, lwd=2, col="grey10")
    mtext("HIV prevalence, age 15-49y", 3, 0.5, font=2, cex=1.2)
    ##
    ## incidence trend
    if(!is.null(fit3))
      plot_incid(fitaggr, fitincrr, fit3, col=cols)
    else
      plot_incid(fitaggr, fitincrr, col=cols)
    if(!is.null(specres))
      lines(as.numeric(names(incid(specres))), incid(specres), lty=2, lwd=2, col="grey10")
    mtext("HIV incidence rate, age 15-49y", 3, 0.5, font=2, cex=1.2)
    ##
    ## r-trend
    if(!is.null(fit3))
      plot_rvec(fitaggr, fitincrr, fit3, col=cols)
    else
      plot_rvec(fitaggr, fitincrr, col=cols)
    mtext("r(t)", 3, 0.5, font=2, cex=1.2)
    ##
    ## vinfl
    ##
    ## sex incrr
    logincrrsex <- log(sapply(fitincrr$param, "[[", "incrr_sex")[1,])
    dens <- density(logincrrsex)
    densCI <- which(dens$x >= quantile(logincrrsex, 0.025) &
                    dens$x <= quantile(logincrrsex, 0.975))
    if(!is.null(fit3)){
      logincrrsex3 <- log(sapply(fit3$param, "[[", "incrr_sex")[1,])
      dens3 <- density(logincrrsex3)
      dens3CI <- which(dens3$x >= quantile(logincrrsex3, 0.025) &
                       dens3$x <= quantile(logincrrsex3, 0.975))
    }
    plot(dens, xlim=c(-0.1, 0.75), col=cols[2],
         main="F:M incidence rate ratio, posterior density", xlab="log F:M IRR")
    polygon(dens$x[c(min(densCI), densCI, max(densCI))], c(0, dens$y[densCI], 0),
            col=transp(cols[2]), border=NA)
    if(!is.null(fit3)){
      lines(dens3, col=cols[3])
      polygon(dens3$x[c(min(dens3CI), dens3CI, max(dens3CI))], c(0, dens3$y[dens3CI], 0),
              col=transp(cols[3]), border=NA)
    }
    segments(x0=mean(log(sapply(fitincrr$param, "[[", "incrr_sex")[1,])), y0=0, y1=6, col=cols[2], lwd=2)
    if(!is.null(fit3))
      segments(x0=mean(log(sapply(fit3$param, "[[", "incrr_sex")[1,])), y0=0, y1=6, col=cols[3], lwd=2)
    segments(x0=log(tail(fitaggr$fp$incrr_sex, 1)), y0=0, y1=6, col=cols[1], lwd=2)
    lines(seq(-0.1, 0.75, 0.01), dnorm(seq(-0.1, 0.75, 0.01), eppasm:::sexincrr.pr.mean, eppasm:::sexincrr.pr.sd), col="darkblue", lty=2)
    segments(x0=eppasm:::sexincrr.pr.mean, y0=0, y1=2, col="darkblue", lwd=2, lty=2)
    ##
    ## age incrr
    xx <- c(1:2, 4:7)
    logincrrage <- estci(sapply(fitincrr$param, "[[", "logincrr_age"))
    if(!is.null(fit3))
      logincrrage3 <- estci(sapply(fit3$param, "[[", "logincrr_age"))
    ## plot men
    plot(1:7+0.5, 1:7, type="n",  xlim=c(1, 8), ylim=c(-2.5, 2),
         xlab="Age group", ylab="log IRR", main = "Male incidence rate ratios (log)", xaxt="n")
    axis(1, 1:7+0.5, paste0(3:9*5, "-", 3:9*5+4))
    abline(h=0, col="grey")
    points(3.5, 0, pch=4, lwd=2.5, col=1, cex=1.2)
    rect(xx+0.1,
         eppasm:::ageincrr.pr.mean[1:6]-qnorm(0.975)*eppasm:::ageincrr.pr.sd,
         xx+0.9,
         eppasm:::ageincrr.pr.mean[1:6]+qnorm(0.975)*eppasm:::ageincrr.pr.sd,
         col=transp("darkblue", 0.1), border=transp("darkblue", 0.5), lty=2)
    ##
    rect(xx+0.1,  logincrrage[xx,3], xx+0.9, logincrrage[xx,4], col=transp(cols[2]), border=NA)
    if(!is.null(fit3))
      rect(xx+0.1,  logincrrage3[xx,3], xx+0.9, logincrrage3[xx,4], col=transp(cols[3]), border=NA)
    defaultincrr <- log(fitaggr$fp$incrr_age[(xx-1)*5+1,1,fitaggr$fp$ss$PROJ_YEARS])
    segments(xx+0.1, defaultincrr, xx+0.9, col=cols[1], lwd=2)
    segments(xx+0.1, eppasm:::ageincrr.pr.mean[1:6], xx+0.9, col=transp("darkblue", 0.5), lwd=1.5, lty=2)
    segments(xx+0.1, logincrrage[xx,1], xx+0.9, col=cols[2], lwd=2)
    if(!is.null(fit3))
      segments(xx+0.1, logincrrage3[xx,1], xx+0.9, col=cols[3], lwd=2)
    ## plot women
    plot(1:7+0.5, 1:7, type="n",  xlim=c(1, 8), ylim=c(-2.5, 2),
         xlab="Age group", ylab="log IRR", main = "Female incidence rate ratios (log)", xaxt="n")
    axis(1, 1:7+0.5, paste0(3:9*5, "-", 3:9*5+4))
    abline(h=0, col="grey")
    points(3.5, 0, pch=4, lwd=2.5, col=1, cex=1.2)
    rect(xx+0.1,
         eppasm:::ageincrr.pr.mean[7:12]-qnorm(0.975)*eppasm:::ageincrr.pr.sd,
         xx+0.9,
         eppasm:::ageincrr.pr.mean[7:12]+qnorm(0.975)*eppasm:::ageincrr.pr.sd,
         col=transp("darkblue", 0.1), border=transp("darkblue", 0.5), lty=2)
    ##
    offset <- nrow(logincrrage)/2
    rect(xx+0.1,  logincrrage[xx+offset,3], xx+0.9, logincrrage[xx++offset,4], col=transp(cols[2]), border=NA)
    if(!is.null(fit3)){
      offset3 <- nrow(logincrrage3)/2
      rect(xx+0.1,  logincrrage3[xx+offset3,3], xx+0.9, logincrrage3[xx++offset3,4], col=transp(cols[3]), border=NA)
    }
    defaultincrr <- log(fitaggr$fp$incrr_age[(xx-1)*5+1,2,fitaggr$fp$ss$PROJ_YEARS])
    segments(xx+0.1, defaultincrr, xx+0.9, col=cols[1], lwd=2)
    segments(xx+0.1, eppasm:::ageincrr.pr.mean[7:12], xx+0.9, col=transp("darkblue", 0.5), lwd=1.5, lty=2)
    segments(xx+0.1, logincrrage[xx+offset,1], xx+0.9, col=cols[2], lwd=2)
    if(!is.null(fit3))
      segments(xx+0.1, logincrrage3[xx+offset3,1], xx+0.9, col=cols[3], lwd=2)
    ##
    mtext(paste0(icountry, ", ", eppmod, "; posterior distribution"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  }
  ##
  ##  Age specific prevalence compared to survey
  ##
  if(2 %in% pages){
    plot_compare_ageprev(fitaggr, fitincrr, fit3, specres, col=cols)
    mtext(paste0(icountry, ", ", eppmod, "; Age-specific prevalence"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  }
  ##
  ##
  ##  Age prevalence trend by age group
  ##
  if(3 %in% pages){
    par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
    ##
    for(iagegr in 1:3){
      for(isex in 1:2){
        strsex <- c("male", "female")[isex]
        stragegr <- c("15-24", "25-34", "35-49")[iagegr]
        aggr.yprev <- estci(fitaggr$agegr3prev[iagegr,isex,,])
        fit.yprev <- estci(fitincrr$agegr3prev[iagegr,isex,,])
        if(!is.null(fit3))
          fit3.yprev <- estci(fit3$agegr3prev[iagegr,isex,,])
        survdat <- subset(agegr3dat, sex == strsex & agegr3==stragegr)
        ##
        xx <- 1998+seq_len(nrow(fit.yprev))
        plot(xx, fit.yprev[,1], type="n", ylim=c(0, 0.05*ceiling(max(survdat$ci_u, aggr.yprev[,1], fit.yprev[,1])/0.05)),
             ylab="", xlab="", main=paste(strsex, stragegr))
        cred.region(xx, t(aggr.yprev[,3:4]), col=transp(cols[1], 0.4))
        cred.region(xx, t(fit.yprev[,3:4]), col=transp(cols[2], 0.4))
        if(!is.null(fit3))
          cred.region(xx, t(fit3.yprev[,3:4]), col=transp(cols[3], 0.4))
        lines(xx, aggr.yprev[,1], col=cols[1], lwd=2)
        lines(xx, fit.yprev[,1], col=cols[2], lwd=2)
        if(!is.null(fit3))
          lines(xx, fit3.yprev[,1], col=cols[3], lwd=2)
        ##
        if(!is.null(specres)){
          ages <- as.character(c(15, 25, 35)[iagegr] + 0:(c(10, 10, 15)[iagegr]-1))
          specres.prev <- colSums(specres$hivpop[ages,isex,as.character(xx)]) / colSums(specres$totpop[ages,isex,as.character(xx)])
          lines(xx, specres.prev, lty=2, lwd=2, col="grey10")
        }
        ##
        points(survdat$year, survdat$prev, pch=19)
        segments(survdat$year, y0=survdat$ci_l, y1=survdat$ci_u)
      }
    }
    mtext(paste0(icountry, ", ", eppmod, "; Prevalence trend by age group"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  }
  ##
  ##
  return(invisible())
}
