plot_prev2 <- function(fit, ..., ylim=NULL, xlim=c(1980, max(as.integer(dimnames(fit$prev)[[1]]))),
                      col="blue", main="", ylab="prevalence"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(fit$prev[,"upper"]))
  xx <- as.integer(dimnames(fit$prev)[[1]])
  plot(xx, fit$prev[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
    
  for(ii in seq_along(dots))
    cred.region(xx, dots[[ii]]$prev[,c("lower", "upper")], col=transp(col[1+ii], 0.3))
  cred.region(xx, fit$prev[,c("lower", "upper")], col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, dots[[ii]]$prev[,"mean"], col=col[1+ii], lwd=1.5)
  lines(xx, fit$prev[,"mean"], col=col[1], lwd=1.5)
}


plot_incid2 <- function(fit, ..., ylim=NULL, xlim=c(1980, max(as.integer(dimnames(fit$incid)[[1]]))),
                       col="blue", main="", ylab="incidence rate"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(fit$incid[,"upper"]))
  xx <- as.integer(dimnames(fit$incid)[[1]])
  plot(xx, fit$incid[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
    
  for(ii in seq_along(dots))
    cred.region(xx, dots[[ii]]$incid[,c("lower", "upper")], col=transp(col[1+ii], 0.3))
  cred.region(xx, fit$incid[,c("lower", "upper")], col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, dots[[ii]]$incid[,"mean"], col=col[1+ii], lwd=1.5)
  lines(xx, fit$incid[,"mean"], col=col[1], lwd=1.5)
}


plot_log_transmrate <- function(fit, ..., ylim=NULL, xlim=c(1980, max(as.integer(dimnames(fit$transmrate)[[1]]))),
                                col="blue", main="", ylab="log transmission rate"){
  if(is.null(ylim))
    ylim <- c(min(log(fit$transmrate[fit$transmrate[,"lower"] > 0, "lower"]))-0.2,
              max(log(fit$transmrate[,"upper"])) + 0.2)
  xx <- as.integer(dimnames(fit$transmrate)[[1]])
  plot(xx, fit$transmrate[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
    
  for(ii in seq_along(dots))
    cred.region(xx, log(dots[[ii]]$transmrate[,c("lower", "upper")]), col=transp(col[1+ii], 0.3))
  cred.region(xx, log(fit$transmrate[,c("lower", "upper")]), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, log(dots[[ii]]$transmrate[,"mean"]), col=col[1+ii], lwd=1.5)
  lines(xx, log(fit$transmrate[,"mean"]), col=col[1], lwd=1.5)
}


plot_incidsexratio <- function(fit, ..., ylim=NULL, xlim=c(1999, max(as.integer(dimnames(fit$incidsexratio)[[1]]))),
                       col="blue", main="", ylab="F:M incidence ratio"){
  if(is.null(ylim))
    ylim <- c(0, max(2.5, 1.1*max(fit$incidsexratio[,"upper"])))
  xx <- as.integer(dimnames(fit$incidsexratio)[[1]])
  plot(xx, fit$incidsexratio[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
    
  for(ii in seq_along(dots))
    cred.region(xx, dots[[ii]]$incidsexratio[,c("lower", "upper")], col=transp(col[1+ii], 0.3))
  cred.region(xx, fit$incidsexratio[,c("lower", "upper")], col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, dots[[ii]]$incidsexratio[,"mean"], col=col[1+ii], lwd=1.5)
  lines(xx, fit$incidsexratio[,"mean"], col=col[1], lwd=1.5)
}



plot_pregprev <- function(fit, ..., likdat=NULL, ylim=NULL, xlim=c(1988, max(as.integer(dimnames(fit$pregprev)[[1]]))),
                          col="blue", main="", ylab="Preg. prevalence"){
  
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
  
  if(is.null(ylim)){
    maxest <- max(fit$pregprev[,"upper"])
    if(!is.null(likdat) && !is.null(likdat$anclik.dat) && length(likdat$anclik.dat$W.lst))
      maxdata <- max(pnorm(unlist(likdat$anclik.dat$W.lst)))
    else
      maxdata <- 0
    ylim <- c(0, 1.1*max(maxest, min(maxdata, 2*maxest)))
  }
  xx <- as.integer(dimnames(fit$incidsexratio)[[1]])
  
  plot(xx, fit$pregprev[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  ##
  if(!is.null(likdat)){
    if(!is.null(likdat$anclik.dat) && length(likdat$anclik.dat$W.lst)){
      with(likdat$anclik.dat, mapply(function(idx, W) points(idx+1970-1, pnorm(W), col=adjustcolor("grey", 0.5), pch=15), anc.idx.lst, W.lst))
      with(likdat$anclik.dat, mapply(function(idx, W) lines(idx+1970-1, pnorm(W), col=adjustcolor("grey", 0.5)), anc.idx.lst, W.lst))
    }
  }
  ##
  for(ii in seq_along(dots))
    cred.region(xx, dots[[ii]]$pregprev[,c("lower", "upper")], col=transp(col[1+ii], 0.3))
  cred.region(xx, fit$pregprev[,c("lower", "upper")], col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, dots[[ii]]$pregprev[,"mean"], col=col[1+ii], lwd=1.5)
  lines(xx, fit$pregprev[,"mean"], col=col[1], lwd=1.5)
  ##
  if(!is.null(likdat$ancrtcens.dat) && length(likdat$ancrtcens.dat$W.ancrt)){
    with(likdat$ancrtcens.dat,
         segments(year, y0=pnorm(qnorm(prev) - qnorm(0.975)*sqrt(v.ancrt)),
                  y1=pnorm(qnorm(prev) + qnorm(0.975)*sqrt(v.ancrt))))
    with(likdat$ancrtcens.dat, points(year, prev, pch=15))
  }
}


plot_artcov15plus <- function(fit, ..., ylim=NULL, xlim=c(2003, max(as.integer(dimnames(fit$artcov15plus)[[1]]))),
                       col="blue", main="", ylab="ART coverage"){
  if(is.null(ylim))
    ylim <- c(0, 1)
  xx <- as.integer(dimnames(fit$artcov15plus)[[1]])
  plot(xx, fit$artcov15plus[,"mean"], type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  dots <- dots[!sapply(dots, is.null)]
    
  for(ii in seq_along(dots))
    cred.region(xx, dots[[ii]]$artcov15plus[,c("lower", "upper")], col=transp(col[1+ii], 0.3))
  cred.region(xx, fit$artcov15plus[,c("lower", "upper")], col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, dots[[ii]]$artcov15plus[,"mean"], col=col[1+ii], lwd=1.5)
  lines(xx, fit$artcov15plus[,"mean"], col=col[1], lwd=1.5)
}


#' @importFrom binom binom.exact
plot_compare_ageprev2 <- function(fit, fit2=NULL, fit3=NULL, specres=NULL, likdat=NULL, ylim=NULL, col=c("navy", "darkred", "forestgreen"), exact_ci=TRUE){
  if(is.null(ylim)){
    if(!is.null(likdat))
      maxdata <- likdat$hhsage.dat$prev
    else
      maxdata <- 0
    ylim <- c(0, 0.05*ceiling(max(c(fit$ageprevdat$upper, 1.3*maxdata))/0.05))
  }
  ####
  survprev <- fit$ageprevdat
  if(!is.null(likdat)){
    survprev <- merge(likdat$hhsage.dat, fit$ageprevdat, by=c("year", "survyear", "sex", "agegr"), all.x=TRUE)
    if(exact_ci)
      survprev[c("ci_l", "ci_u")] <- with(survprev, binom.exact(x_eff, n_eff))[c("lower", "upper")]
  }
  survprev$survyear <- with(survprev, factor(survyear, levels(survyear)[order(as.integer(substr(levels(survyear), 1, 4)))]))
  survprev <- split(survprev, factor(survprev$survyear))  
  ##
  if(!is.null(fit2))
    survprev2 <- split(fit2$ageprevdat, factor(fit2$ageprevdat$survyear))
  if(!is.null(fit3))
    survprev3 <- split(fit3$ageprevdat, factor(fit3$ageprevdat$survyear))
  ##
  par(mfrow=c(4,2), mar=c(2, 3, 2, 1), tcl=-0.25, mgp=c(2, 0.5, 0), las=1, cex=1)
  for(isurv in names(survprev))
    for(isex in c("male", "female")){
      sp <- subset(survprev[[isurv]], sex==isex & as.integer(agegr) %in% 3:11)
      if(!is.null(fit2))
        sp2 <- subset(survprev2[[isurv]], sex==isex & as.integer(agegr) %in% 3:11)
      if(!is.null(fit3))
        sp3 <- subset(survprev3[[isurv]], sex==isex & as.integer(agegr) %in% 3:11)
      ##
      xx <- as.integer(sp$agegr)
      main <- if(!is.null(sp$eppregion))
                paste0(sp$country[1], " ", gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", sp$eppregion[1], perl=TRUE), " ", survprev[[isurv]]$survyear[1], ", ", isex)
              else
                paste0(sp$country[1], " ", survprev[[isurv]]$survyear[1], ", ", isex)
      plot(xx+0.5, sp$prev, type="n", xlim=c(4, 12), ylim=ylim, xaxt="n",
           main=main, xlab="", ylab="")
      axis(1, xx+0.5, sp$agegr)
      ##
      rect(xx+0.05, sp$lower, xx+0.95, sp$upper,
           col=transp(col[1]), border=NA)
      segments(xx+0.05, sp$mean, xx+0.95, col=col[1], lwd=2)
      ##
      if(!is.null(fit2)){
        rect(xx+0.05, sp2$lower, xx+0.95, sp2$upper,
             col=transp(col[2]), border=NA)
        segments(xx+0.05, sp2$mean, xx+0.95, col=col[2], lwd=2)
      }
      if(!is.null(fit3)){
        rect(xx+0.05, sp3$lower, xx+0.95, sp3$upper,
             col=transp(col[3]), border=NA)
        segments(xx+0.05, sp3$mean, xx+0.95, col=col[3], lwd=2)
      }
      ##
      if(!is.null(specres)){
        csex <- sub("(\\b[a-z]{1})", "\\U\\1" , isex, perl=TRUE)
        stryear <- as.character(survprev[[isurv]]$year[1])
        specres.prev <- tapply(specres$hivpop[as.character(15:54), csex, stryear], rep(3:10, each=5), sum) / tapply(specres$totpop[as.character(15:54), csex, stryear], rep(3:10, each=5), sum)
        segments(4:11+0.1, specres.prev, 4:11+0.9, lty=3, col="grey10", lwd=2)
      }
      if(exists("prev", sp)){
        points(xx+0.5, sp$prev, pch=19)
        segments(x0=xx+0.5, y0=sp$ci_l, y1=sp$ci_u)
      }
    }
  ##
  return(invisible())
}
