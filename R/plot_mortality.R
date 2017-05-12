plot_q4515 <- function(icountry, fit=NULL, fit2=NULL, specres=NULL, plotdata=TRUE, plotgam=TRUE){
  ## device: quartz(h=4, w=8.5, pointsize=12, type="cairo")
  par(mfrow=c(1, 2), las=1, mgp=c(2.0, 0.5, 0), tcl=-0.25, mar=c(2.0, 3.0, 2, 0.5), cex=0.8)
  ##
  for(isex in c("female", "male")){
    dat <- subset(sibq4515, sex==isex & country == icountry)
    cols <- RColorBrewer::brewer.pal(3, "Set1")
    plot(dat$period, dat$q4515, type="n",
         ylim=c(0, 0.9), xlim=c(1985, 2015), ylab=expression(""[45]~q[15]), xlab="", main=paste(icountry, c("male"="Men", "female"="Women")[isex]))
    if(plotdata)
      with(subset(adjq4515, sex == isex & country == icountry), points(period, q4515, pch=23, col=transp("grey", 0.4), bg=transp("grey", 0.4)))
    if(plotgam){
      cred.region(dat$period, t(dat[c("ci_l", "ci_u")]), col=transp(cols[2], 0.2))
      lines(dat$period, dat$q4515, lwd=2, lty=1, col=cols[2])
    }
    ##
    if(!is.null(specres)){
      lines(1983:2015, calc_nqx(specres, 45, 15)[c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=2, lwd=3)
      lines(1983:2015, calc_nqx(specres, 45, 15, nonhiv=TRUE)[c("male"="Male", "female"="Female")[isex], as.character(1983:2015)], lty=3, lwd=3, col="grey30")
    }
    
    if(!is.null(fit)){
      est <- estci(fit$q4515[isex, as.character(1983:2015),])
      cred.region(as.integer(row.names(est)), t(est[,c("lower", "upper")]), col=transp("darkred", 0.4))
      lines(1983:2015, est[,"mean"], col="darkred", lwd=2)
      ## nonhiv
      est <- estci(fit$natq4515[isex, as.character(1983:2015),])
      cred.region(as.integer(row.names(est)), t(est[,c("lower", "upper")]), col=transp("darkred", 0.4))
      lines(1983:2015, est[,"mean"], col="darkred", lwd=1, lty=2)
    }
    if(!is.null(fit2)){
      est <- estci(fit2$q4515[isex, as.character(1983:2015),])
      cred.region(as.integer(row.names(est)), t(est[,c("lower", "upper")]), col=transp("darkolivegreen", 0.3))
      lines(1983:2015, est[,"mean"], col="darkolivegreen", lwd=2)
      ## nonhiv
      est <- estci(fit2$natq4515[isex, as.character(1983:2015),])
      cred.region(as.integer(row.names(est)), t(est[,c("lower", "upper")]), col=transp("darkolivegreen", 0.4))
      lines(1983:2015, est[,"mean"], col="darkolivegreen", lwd=1, lty=2)
    }
    
    ## if(plothhd){
    ##   with(subset(hhd_q4515, country == icountry & sex == isex), segments(year, y0=ci_l, y1=ci_u, lwd=2))
    ##   with(subset(hhd_q4515, country == icountry & sex == isex), points(year, q4515, pch=19))
    ## }
    ## if(plotaggr & !plotagemod)
    ##   legend("topleft", c("Sibling history", "Default IRRs"), pch=15, pt.cex=2.2, col=transp(c(cols[2], "darkolivegreen")), bty="n")
    ## if(plotagemod)
    ##   legend("topleft", c("Sibling history", "Default IRRs", "Estimated IRRs", "UNAIDS 2016"), pch=c(15, 15, 15, NA), pt.cex=2.2, col=c(transp(c(cols[2], "darkolivegreen", "darkred")), 1), lty=c(NA, NA, NA, 2), lwd=c(NA, NA, NA, 2), bty="n")
    ## legend("topleft", c("Sibling history", "Estimated IRRs", "UNAIDS 2016", "UNAIDS '16 non-HIV"), pch=c(15, 15, NA, NA), pt.cex=2.2, col=c(transp(c(cols[2], "darkred")), 1, "grey30"), lty=c(NA, NA, 2, 3), lwd=c(NA, NA, 2, 2), bty="n", cex=0.85)
  }
  return(invisible())
}


plot_agemx <- function(icountry, fit=NULL, fit2=NULL, specres=NULL, plotgam=TRUE,
                       cols=RColorBrewer::brewer.pal(3, "Set1")[c(2,1,3)], periods=NULL, smoothed=FALSE){
  ##
  ## device:: quartz(h=4, w=7.5, pointsize=12)
  ##
  par(mfrow=c(2, 4), las=1, mgp=c(1.5, 0.5, 0), tcl=-0.25, oma=c(0, 3.0, 2.2, 0.5), mar=c(0.5, 0.5, 0.25, 0.25), cex=0.8)
  layout(rbind(1:4, 5:8, 0, 9), h=c(1,1,0.1, 0.2))
  ##
  sibest <- subset(sibagemx, country == icountry)
  ##
  if(is.null(periods))
    periods <- round(seq(max(min(sibest$period)+2, 1988), min(max(sibest$period), 2015), length.out=4))
                                                      
  
  ## 
  for(isex in c("female", "male")){

    for(per in periods){
      dat <- subset(sibagemx, sex == isex & period == per & agegr %in% 15:59 & country == icountry)
      with(dat, plot(agegr, lmx, type="n", ylim=c(-7, -2.5), ylab="", xlab="", main="", xaxt="n", yaxt="n"))
      mtext(paste(c(male="Men", female="Women")[isex], per), 3, -1.1, adj=0.05, font=2, cex=0.8)
      yaxl <- c(1, 2, 5, 10, 25, 50)
      yax <- log(yaxl/1000)
      axis(1, labels=FALSE); axis(2, yax, labels=FALSE)
      if(isex == "male"){axis(1, tick=FALSE); mtext("Age", 1, 1.5, cex=0.8)}
      if(per == periods[1]){axis(2, yax, yaxl, tick=FALSE)}
      if(plotgam){
        with(dat, cred.region(agegr, c(lmx) + (qnorm(0.975)*se.lmx %*% t(c(-1,1))), col=transp(cols[1], 0.2)))
        with(dat, lines(agegr, lmx, col=cols[1], lwd=2, lty=1))
      }

      if(!is.null(specres)){
        lines(15:59, log(agemx(specres)[as.character(15:59), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=2, lwd=2)
        lines(15:59, log(agemx(specres, nonhiv=TRUE)[as.character(15:59), c("male"="Male", "female"="Female")[isex], as.character(per)]), lty=3, lwd=2, col="grey30")
      }

      if(!is.null(fit)){
        logmx <- estci(log(fit$agemx[, isex, as.character(per),]))
        if(smoothed)
          est <- apply(logmx, 2, function(y) predict(smooth.spline(as.integer(row.names(logmx)), y, spar=0.5), 15:59)$y)
        else
          est <- logmx[as.character(15:59),]
        cred.region(15:59, est[,c("lower", "upper")], col=transp(cols[2], 0.3))
        lines(15:59, est[,"mean"], col=cols[2], lwd=2)

        ## non-HIV mortality
        lognatmx <- estci(log(fit$natagemx[, isex, as.character(per),]))
        est <- lognatmx[as.character(15:59),]
        cred.region(15:59, est[,c("lower", "upper")], col=transp(cols[2], 0.3))
        lines(15:59, est[,"mean"], col=cols[2], lty=3, lwd=1.5)
      }

      if(!is.null(fit2)){
        logmx <- estci(log(fit2$agemx[, isex, as.character(per),]))
        if(smoothed)
          est <- apply(logmx, 2, function(y) predict(smooth.spline(as.integer(row.names(logmx)), y, spar=0.5), 15:59)$y)
        else
          est <- logmx[as.character(15:59),]
        cred.region(15:59, est[,c("lower", "upper")], col=transp(cols[3], 0.3))
        lines(15:59, est[,"mean"], col=cols[3], lwd=2)

        ## non-HIV mortality
        lognatmx <- estci(log(fit2$natagemx[, isex, as.character(per),]))
        est <- lognatmx[as.character(15:59),]
        cred.region(15:59, est[,c("lower", "upper")], col=transp(cols[3], 0.3))
        lines(15:59, est[,"mean"], col=cols[3], lty=3, lwd=1.5)
      }
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
  ## ##
  return(invisible())
}
