

plot_output <- function(fit, fit2=NULL, fit3=NULL, specres=NULL, cols=c("navy", "darkred", "darkgreen"), pages=1:3, likdat=NULL,
                        title="",
                        agegr3dat=NULL, age15to49dat=NULL){
  if(1 %in% pages){
    graphics::par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
    ##
    ## prevalence trend
    plot_prev2(fit, fit2, fit3, col=cols)
    if(!is.null(specres))
      graphics::lines(as.numeric(names(prev(specres))), prev(specres), lty=2, lwd=2, col="grey10")
    if(!is.null(age15to49dat)){
      if(!is.null(age15to49dat$prev_spec17)){
          graphics::points(age15to49dat$year, age15to49dat$prev, pch=15, col="grey40")
          graphics::points(age15to49dat$year, age15to49dat$prev_spec17, pch=19)
          graphics::segments(age15to49dat$year, y0=age15to49dat$ci_l_spec17, y1=age15to49dat$ci_u_spec17)
        } else {
          graphics::points(age15to49dat$year, age15to49dat$prev, pch=19)
          graphics::segments(age15to49dat$year, y0=age15to49dat$ci_l, y1=age15to49dat$ci_u)
        }
    }
    graphics::mtext("HIV prevalence, age 15-49y", 3, 0.5, font=2, cex=1.1)
    ##
    ## incidence trend
    plot_incid2(fit, fit2, fit3, col=cols)
    if(!is.null(specres))
      graphics::lines(as.numeric(names(incid(specres))), incid(specres), lty=2, lwd=2, col="grey10")
    graphics::mtext("HIV incidence rate, age 15-49y", 3, 0.5, font=2, cex=1.1)
    ##
    ## tranmission rate
    plot_log_transmrate(fit, fit2, fit3, col=cols)
    if(!is.null(specres))
      graphics::lines(as.numeric(names(incid(specres))), log(incid(specres) / prev(specres)), lty=2, lwd=2, col="grey10")
    graphics::mtext("Log transmission rate", 3, 0.5, font=2, cex=1.1)
    ## 
    ## Sex ratio of incidence
    plot_incidsexratio(fit, fit2, fit3, col=cols)
    if(!is.null(specres))
      graphics::lines(as.numeric(names(incid_sexratio(specres))), incid_sexratio(specres), lty=2, lwd=2, col="grey10")
    graphics::mtext("Incidence sex ratio (F:M, 15-49y)", 3, 0.5, font=2, cex=1.1)
    ##
    ## ART coverage age 15+
    plot_artcov15plus(fit, fit2, fit3, col=cols)
    if(!is.null(specres) && exists("artnum.m", specres))
      graphics::lines(as.numeric(names(artcov15plus(specres))), artcov15plus(specres), lty=2, lwd=2, col="grey10")
    graphics::mtext("ART coverage, age 15+", 3, 0.5, font=2, cex=1.1)
    graphics::mtext(paste0(title, "; posterior distribution"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
    ##
    ## Pregnant women prevalence
    if(!is.null(fit$pregprev)){
      plot_pregprev(fit, fit2, fit3, col=cols, likdat=likdat)
      if(!is.null(specres))
        graphics::lines(as.numeric(names(fnPregPrev(specres))), fnPregPrev(specres), lty=2, lwd=2, col="grey10")
      graphics::mtext("Prevalence among pregnant women", 3, 0.5, font=2, cex=1.1)
    }
  }
  ##
  ##  Age specific prevalence compared to survey
  ##
  if(2 %in% pages){
    plot_compare_ageprev2(fit, fit2, fit3, specres, col=cols, likdat=likdat)
    graphics::mtext(paste0(title, "; Age-specific prevalence"), 3, 0.5, outer=TRUE, font=2, cex=1.3)
  }
  ##
  ##
  ##  Age prevalence trend by age group
  ##
  if(3 %in% pages){
    graphics::par(oma=c(1, 1, 2.5, 1), mfrow=c(3, 2), mar=c(3.5, 3.5, 3, 1), tcl=-0.25, mgp=c(2.5, 0.5, 0), cex=1, las=1)
    ##
    for(iagegr in 1:3){
      for(isex in 1:2){
        strsex <- c("Male", "Female")[isex]
        stragegr <- c("15-24", "25-34", "35-49")[iagegr]
        fit.yprev <- fit$agegr3prev[stragegr,strsex,,]
        if(!is.null(fit2))
          fit2.yprev <- fit2$agegr3prev[stragegr,strsex,,]
        else
          fit2.yprev <- NULL
        if(!is.null(fit3))
          fit3.yprev <- fit3$agegr3prev[stragegr,strsex,,]
        else
          fit3.yprev <- NULL
        if(!is.null(agegr3dat))
          survdat <- subset(agegr3dat, sex == tolower(strsex) & agegr3==stragegr)
        else
          survdat <- NULL
        ##
        xx <- as.integer(rownames(fit.yprev))
        maxprev <- max(survdat$ci_u, fit.yprev[,"upper"], fit2.yprev[,"upper"], fit3.yprev[,"upper"], na.rm=TRUE)
        plot(xx, fit.yprev[,"mean"], type="n",
             ylim=c(0, 0.05*ceiling(maxprev/0.05)),
             ylab="", xlab="", main=paste(strsex, stragegr))
        cred.region(xx, t(fit.yprev[,c("lower", "upper")]), col=transp(cols[1], 0.4))
        if(!is.null(fit2))
          cred.region(xx, t(fit2.yprev[,c("lower", "upper")]), col=transp(cols[2], 0.4))
        if(!is.null(fit3))
          cred.region(xx, t(fit3.yprev[,c("lower", "upper")]), col=transp(cols[3], 0.4))
        graphics::lines(xx, fit.yprev[,"mean"], col=cols[1], lwd=2)
        if(!is.null(fit2))
          graphics::lines(xx, fit2.yprev[,"mean"], col=cols[2], lwd=2)
        if(!is.null(fit3))
          graphics::lines(xx, fit3.yprev[,"mean"], col=cols[3], lwd=2)
        ##
        if(!is.null(specres)){
          ages <- as.character(c(15, 25, 35)[iagegr] + 0:(c(10, 10, 15)[iagegr]-1))
          specres.prev <- colSums(specres$hivpop[ages,isex,as.character(xx)]) / colSums(specres$totpop[ages,isex,as.character(xx)])
          graphics::lines(xx, specres.prev, lty=2, lwd=2, col="grey10")
        }
        ##
        if(!is.null(survdat)){
          if(!is.null(survdat$prev_spec17)){
            graphics::points(survdat$year, survdat$prev, pch=15, col="grey40")
            graphics::points(survdat$year, survdat$prev_spec17, pch=19)
            graphics::segments(survdat$year, y0=survdat$ci_l_spec17, y1=survdat$ci_u_spec17)
          } else {
            graphics::points(survdat$year, survdat$prev, pch=19)
            graphics::segments(survdat$year, y0=survdat$ci_l, y1=survdat$ci_u)
          }
        }
      }
    }
    graphics::mtext(paste0(title, "\nPrevalence trend by age group"), 3, 0, outer=TRUE, font=2, cex=1.3)
  }
  ##
  ##
  return(invisible())
}
