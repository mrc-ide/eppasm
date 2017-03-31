cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), border=NA, ...)

transp <- function(col, alpha=0.5)
  return(replace(apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha)), is.na(col), NA))

estci <- function(x){val <- cbind(rowMeans(x), t(apply(x, 1, quantile, c(0.5, 0.025, 0.975)))); colnames(val) <- c("mean", "median", "lower", "upper"); val}

plot_compare_ageprev <- function(fit, fit2=NULL, fit3=NULL, specres=NULL, ylim=NULL, col=c("grey30", "darkred", "forestgreen")){
  if(is.null(ylim))
    ylim <- c(0, 0.05*ceiling(max(fit$likdat$hhsage.dat$ci_u)/0.05))
  ####
  survprev <- data.frame(fit$likdat$hhsage.dat,
                         estci(fit$ageprevdat))
  survprev <- split(survprev, factor(survprev$year))
  ##
  if(!is.null(fit2)){
    survprev2 <- data.frame(fit2$likdat$hhsage.dat,
                            estci(fit2$ageprevdat))
    survprev2 <- split(survprev2, factor(survprev2$year))
  }
  if(!is.null(fit3)){
    survprev3 <- data.frame(fit3$likdat$hhsage.dat,
                            estci(fit3$ageprevdat))
    survprev3 <- split(survprev3, factor(survprev3$year))
  }
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
      main <- if(!is.null(sp$region))
                paste0(sp$country[1], " ", gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", sp$region[1], perl=TRUE), " ", survprev[[isurv]]$survyear[1], ", ", isex)
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
        specres.prev <- tapply(specres$hivpop[as.character(15:54), csex, isurv], rep(3:10, each=5), sum) / tapply(specres$totpop[as.character(15:54), csex, isurv], rep(3:10, each=5), sum)
        segments(4:11+0.1, specres.prev, 4:11+0.9, lty=3, col="grey10", lwd=2)
      }
      points(xx+0.5, sp$prev, pch=19)
      segments(x0=xx+0.5, y0=sp$ci_l, y1=sp$ci_u)
    }
  ####
  return(invisible())
}



plot_prev <- function(fit, ..., ylim=NULL, xlim=c(1980, 2016), col="blue", main="", ylab="prevalence"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(apply(fit$prev, 1, quantile, 0.975)))
  xx <- fit$fp$ss$proj_start-1+1:fit$fp$ss$PROJ_YEARS
  plot(xx, rowMeans(fit$prev), type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  for(ii in seq_along(dots))
    cred.region(xx, apply(dots[[ii]]$prev, 1, quantile, c(0.025, 0.975)), col=transp(col[1+ii], 0.3))
  cred.region(xx, apply(fit$prev, 1, quantile, c(0.025, 0.975)), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, rowMeans(dots[[ii]]$prev), col=col[1+ii], lwd=1.5)
  lines(xx, rowMeans(fit$prev), col=col[1], lwd=1.5)
  ##
  points(fit$likdat$hhslik.dat$year, fit$likdat$hhslik.dat$prev, pch=20)
  segments(fit$likdat$hhslik.dat$year,
           y0=pnorm(fit$likdat$hhslik.dat$W.hhs - qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs),
           y1=pnorm(fit$likdat$hhslik.dat$W.hhs + qnorm(0.975)*fit$likdat$hhslik.dat$sd.W.hhs))
}

plot_incid <- function(fit, ..., ylim=NULL, xlim=c(1980, 2016), col="blue", main="", ylab="incidence rate"){
  if(is.null(ylim))
    ylim <- c(0, 1.1*max(apply(fit$incid, 1, quantile, 0.975)))
  xx <- fit$fp$ss$proj_start-1+1:fit$fp$ss$PROJ_YEARS
  plot(xx, rowMeans(fit$incid), type="n", ylim=ylim, xlim=xlim,
       ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  dots <- list(...)
  for(ii in seq_along(dots))
    cred.region(xx, apply(dots[[ii]]$incid, 1, quantile, c(0.025, 0.975)), col=transp(col[1+ii], 0.3))
  cred.region(xx, apply(fit$incid, 1, quantile, c(0.025, 0.975)), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, rowMeans(dots[[ii]]$incid), col=col[1+ii],lwd=1.5)
  lines(xx, rowMeans(fit$incid), col=col, lwd=1.5)
  ##
  if(exists("hhsincid.dat", where=fit$likdat))
    with(fit$likdat$hhsincid.dat,{ points(year, incid, pch=20);
      segments(year, y0=exp(log_incid-qnorm(0.975)*log_incid.se), y1=exp(log_incid+qnorm(0.975)*log_incid.se))})
}

plot_rvec <- function(fit, ..., ylim=NULL, xlim=c(1980, 2016), col="blue"){
  dots <- list(...)
  fit$rvec <- sapply(seq_len(ncol(fit$rvec)), function(i) replace(fit$rvec[,i], fit$fp$proj.steps < fit$param[[i]]$tsEpidemicStart, NA))
  idx <- seq_len(length(fit$fp$proj.steps)-1)
  xx <- fit$fp$proj.steps[idx]
  for(ii in seq_along(dots))
    dots[[ii]]$rvec <- sapply(seq_len(ncol(dots[[ii]]$rvec)), function(i) replace(dots[[ii]]$rvec[,i], dots[[ii]]$fp$proj.steps < dots[[ii]]$param[[i]]$tsEpidemicStart, NA))
  if(is.null(ylim))
    ylim <- c(0, min(quantile(apply(fit$rvec, 1, quantile, 0.975, na.rm=TRUE), 0.975, na.rm=TRUE), 5))
  plot(xx, rowMeans(fit$rvec, na.rm=TRUE)[idx], type="n", ylim=ylim, xlim=xlim, ylab="r(t)", yaxt="n", xlab="")
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  for(ii in seq_along(dots))
    cred.region(xx, apply(dots[[ii]]$rvec[idx,], 1, quantile, c(0.025, 0.975), na.rm=TRUE), col=transp(col[1+ii], 0.3))
  cred.region(xx, apply(fit$rvec[idx,], 1, quantile, c(0.025, 0.975), na.rm=TRUE), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, rowMeans(dots[[ii]]$rvec[idx,], na.rm=TRUE)[idx], col=col[1+ii])
  lines(xx, rowMeans(fit$rvec[idx,], na.rm=TRUE)[idx], col=col[1])
  return(invisible())
}

plot_ancprev <- function(fit, ..., data=fit$likdat, ylim=NULL, xlim=c(1980, 2016), col="blue", main="", ylab="ANC prevalence"){
  
  ## Calculate parameters for bias term
  param <- apply(fit$resample, 1, fnCreateParam, fit$fp)
  ancrtcens.bias <- sapply(param, "[[", "ancrtcens.bias")
  ancrtcens.vinfl <- sapply(param, "[[", "ancrtcens.vinfl")
  ancrt.prev <- pnorm(sweep(qnorm(fit$pregprev), 2, ancrtcens.bias, "+"))

  dots <- list(...)
  dots_param <- lapply(dots, function(fit) apply(fit$resample, 1, fnCreateParam, fit$fp))
  dots_ancrtcens.bias <- lapply(dots_param, sapply, "[[", "ancrtcens.bias")
  dots_ancrtcens.vinfl <- lapply(dots_param, sapply, "[[", "ancrtcens.vinfl")
  dots_ancrt.prev <- mapply(function(fit, bias) pnorm(sweep(qnorm(fit$pregprev), 2, bias, "+")), dots, dots_ancrtcens.bias, SIMPLIFY=FALSE)

  if(is.null(ylim))
    ylim <- c(0, 1.1*max(apply(ancrt.prev, 1, quantile, 0.975)))
  xx <- fit$fp$ss$proj_start-1+1:fit$fp$ss$PROJ_YEARS
  plot(xx, rowMeans(fit$pregprev), type="n", ylim=ylim, xlim=xlim, ylab=ylab, xlab="", yaxt="n", xaxt="n", main=main)
  axis(1, labels=TRUE)
  axis(2, labels=TRUE)
  for(ii in seq_along(dots))
    cred.region(xx, apply(dots_ancrt.prev[[ii]], 1, quantile, c(0.025, 0.975)), col=transp(col[1+ii], 0.3))
  cred.region(xx, apply(ancrt.prev, 1, quantile, c(0.025, 0.975)), col=transp(col[1], 0.3))
  for(ii in seq_along(dots))
    lines(xx, rowMeans(dots_ancrt.prev[[ii]]), col=col[1+ii], lwd=1.5)
  lines(xx, rowMeans(ancrt.prev), col=col[1], lwd=1.5)
  ##
  if(!is.null(data$ancrtcens.dat)){
    vinfl <- if(!is.null(fit$fp$ancrtcens.vinfl)) fit$fp$ancrtcens.vinfl else mean(ancrtcens.vinfl)
    with(data$ancrtcens.dat, segments(x0=year, col="grey50",
                                      y0=pnorm(qnorm(prev) - qnorm(0.975)*sqrt(v.ancrt+vinfl)),
                                      y1=pnorm(qnorm(prev) + qnorm(0.975)*sqrt(v.ancrt+vinfl))))
    with(data$ancrtcens.dat, segments(year,
                                      y0=pnorm(qnorm(prev) - qnorm(0.975)*sqrt(v.ancrt)),
                                      y1=pnorm(qnorm(prev) + qnorm(0.975)*sqrt(v.ancrt))))
    with(data$ancrtcens.dat, points(year, prev, pch=20))
  }
}
