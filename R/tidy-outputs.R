tidy_output <- function(fit, modlab, country=NA, eppregion=NA){
  print(paste(country, eppregion))

  ss <- fit$fp$ss
  
  ## pregprev predictions
  pred <- rbind(data.frame(year = 1980:2018, age = 15, agspan = 35),
                expand.grid(year = 1990:2018, age = 3:9*5, agspan =5))
  pred$agegr <- with(pred, paste0(age, "-", age+agspan-1))
  pred$agegr <- relevel(factor(pred$agegr), "15-49")
  pred$aidx <- pred$age - ss$AGE_START + 1
  pred$yidx <- pred$year - ss$proj_start+ 1

  ## agegr3 predictions
  agegr3 <- merge(expand.grid(year = 1999:2018, sex=c("male", "female"),
                              age = c(15, 25, 35)),
                  data.frame(age = c(15, 25, 35), agspan = c(10, 10, 15)))
  agegr3$agegr3 <- with(agegr3, paste0(age, "-", age+agspan-1))
  agegr3$sidx <- match(agegr3$sex, c("both", "male", "female")) - 1L
  agegr3$aidx <- agegr3$age - ss$AGE_START + 1
  agegr3$yidx <- agegr3$year - ss$proj_start + 1

  fit <- simfit(fit, pregprev = pred, ageprevdat = fit$likdat$hhs.dat,
                agegr3=agegr3, aidsdeaths = TRUE,
                artcov = TRUE, ancartcov=TRUE)

  ## Assemble results

  year <- fit$fp$ss$proj_start + 1:fit$fp$ss$PROJ_YEARS - 1L

  fit$artcov[is.na(fit$artcov)] <- 0
  fit$ancartcov[is.na(fit$ancartcov)] <- 0

  core <- rbind(data.frame(country = country,
                           eppregion = eppregion,
                           year = year,
                           indicator="prev",
                           model = modlab,
                           estci2(fit$prev)),
                data.frame(country = country,
                           eppregion = eppregion,
                           year = year,
                           indicator = "incid",
                           model = modlab,
                           estci2(fit$incid)),
                data.frame(country = country,
                           eppregion = eppregion,
                           year = year,
                           indicator = "aidsdeaths15pl",
                           model = modlab,
                           estci2(fit$aidsdeaths)),
                subset(data.frame(country = country,
                                  eppregion = eppregion,
                                  year = year,
                                  indicator = "artcov15pl",
                                  model = modlab,
                                  estci2(fit$artcov)),
                       year >= 2000),
                subset(data.frame(country = country,
                                  eppregion = eppregion,
                                  year = year,
                                  indicator = "pregartcov",
                                  model = modlab,
                                  estci2(fit$ancartcov)),
                       year >= 2000))

  rvecidx <- seq_len(ss$PROJ_YEARS-1)*ss$hiv_steps_per_year

  core <- rbind(core,
                data.frame(country = country,
                           eppregion = eppregion,
                           year = year[-1],
                           indicator = "r(t)",
                           model = modlab,
                           estci2(fit$rvec[rvecidx, ])),
                data.frame(country = country,
                           eppregion = eppregion,
                           year = year[-1],
                           indicator = "log r(t)",
                           model = modlab,
                           estci2(log(fit$rvec[rvecidx, ]))))

  core <- rbind(core,
                data.frame(country = country,
                           eppregion = eppregion,
                           year = dimnames(fit$agegr3incid)[[3]],
                           model = modlab,
                           indicator = "incidence sex ratio",
                           estci2(fit$agegr3incid["15-49", "female",,] /
                                  fit$agegr3incid["15-49", "male",,])))

  pregprev <- data.frame(country = country,
                         eppregion = eppregion,
                         pred, model=modlab, estci2(fit$pregprev))

  agegr3 <- data.frame(country = country,
                       eppregion = eppregion,
                       agegr3, model = modlab, estci2(fit$agegr3prev))


  ageincid <- reshape2::melt(estci2(fit$agegr3incid),
                        varnames = c("agegr", "sex", "year", "outcome"))
  ageincid <- data.frame(country = country,
                         eppregion = eppregion,
                         indicator = "incidence",
                         model=modlab,
                         reshape2::dcast(ageincid, ... ~ outcome, value.var="value"))

  ageinfections <- reshape2::melt(estci2(fit$agegr3incid),
                                  varnames = c("agegr", "sex", "year", "outcome"))
  ageinfections <- data.frame(country = country,
                              eppregion = eppregion,
                              indicator = "infections",
                              model=modlab,
                              reshape2::dcast(ageinfections, ... ~ outcome, value.var="value"))

  relincid <- reshape2::melt(fit$relincid,
                             varnames = c("agegr", "sex", "year", "outcome"))
  relincid <- data.frame(country = country,
                         eppregion = eppregion,
                         indicator = "relincid",
                         model=modlab,
                         reshape2::dcast(relincid, ... ~ outcome, value.var="value"))

  ageprevdat <- data.frame(country = country,
                           eppregion = eppregion,
                           model = modlab,
                           fit$likdat$hhs.dat[c("year", "sex", "agegr", "n", "prev", "se", "ci_l", "ci_u")],
                           estci2(fit$ageprevdat))
  
  out <- list(country = country,
              eppregion = eppregion,
              core = core,
              ageprevdat=ageprevdat,
              pregprev = pregprev,
              agegr3 = agegr3,
              ageincid = ageincid,
              ageinfections = ageinfections,
              relincid = relincid)
}


combine_tidy <- function(...){

  dots <- list(...)
  vars <- setdiff(names(dots[[1]]), c("country", "eppregion"))

  out <- dots[[1]][c("country", "eppregion")]
  c(out, Map(rbind, ...))
}


plot_output <- function(out, data_color = "grey20", th = theme_light()){

  ## Core indicators
  out$pregprev$indicator <- "Prevalence among pregnant women"
  df <- rbind(subset(out$core, year %in% 1990:2018),
              subset(out$pregprev, agegr == "15-49",
                     c(country, eppregion, year, indicator, model, mean, se, median, lower, upper)))

  obs <- data.frame(indicator = "Prevalence 15-49y",
                    subset(prev_15to49_eppregion,
                           country == out$country &
                           eppregion == out$eppregion & used))

  ggCore <- ggplot(df) +
    geom_point(aes(year, prev), data=obs, col=data_color) +
    geom_errorbar(aes(year, ymin=ci_l, ymax=ci_u), data=obs, width=0, color=data_color) +
    geom_line(aes(year, mean, color=model)) +
    geom_ribbon(aes(year, ymin=lower, ymax=upper, fill=model), alpha=0.3) +
    facet_wrap(~indicator, ncol=2, scales="free") +
    scale_y_continuous(element_blank()) +
    th +
    theme(legend.position = "bottom",
          plot.margin=unit(c(1, 0.5, 1, 0.5), "in")) +
    ggtitle(paste0(out$country, " - ", out$eppregion))


  ## Broad age group prevalence trends

  df <- out$agegr3
  obs <- subset(prev_agegr3sex_eppregion,
                country == out$country &
                eppregion == out$eppregion &
                agegr3 %in% out$agegr3$agegr3 &
                !(country == "Zambia" & survyear %in% c("2013-14", "2013-14retest")))

  ggAgegr3 <- ggplot(df) +
    geom_point(aes(year, prev), data=obs, col=data_color) +
    geom_errorbar(aes(year, ymin=ci_l, ymax=ci_u), data=obs, width=0, color=data_color) +
    geom_line(aes(year, mean, color=model)) +
    geom_ribbon(aes(year, ymin=lower, ymax=upper, fill=model), alpha=0.3) +
    facet_grid(agegr3 ~ sex) +
    th +
    theme(legend.position = "bottom",
          plot.margin=unit(c(1, 0.5, 1, 0.5), "in")) +
    ggtitle(paste0(out$country, " - ", out$eppregion), "Prevalence trend by broad age groups")

  ## ANC prevalence

  df <- out$pregprev
  df$agegr <- relevel(factor(df$agegr), "15-49")

  obs <- subset(ancrtcens, country == out$country & eppregion == out$eppregion)
  if(nrow(obs))
    obs$agegr <- with(obs, paste0(age, "-", age+agspan-1))
  if(nrow(obs))
    obs$agegr <- factor(obs$agegr, levels(df$agegr))
  x.ancrt <- (obs$prev*obs$n+0.5)/(obs$n+1)
  obs$W.ancrt <- qnorm(x.ancrt)
  obs$v.ancrt <- 2*pi*exp(obs$W.ancrt^2)*x.ancrt*(1-x.ancrt)/obs$n
  obs$ci_l <- with(obs, pnorm(W.ancrt - qnorm(0.975) * sqrt(v.ancrt)))
  obs$ci_u <- with(obs, pnorm(W.ancrt + qnorm(0.975) * sqrt(v.ancrt)))

  ancsiteobs <- subset(ancsitedata, country == out$country &
                                    eppregion == out$eppregion &
                                    used)
  ancsiteobs$agegr <- factor(ancsiteobs$agegr, levels(df$agegr))
  
  ggPregPrev <- ggplot(df) +
    geom_point(aes(year, prev), data=obs, col=data_color) +
    geom_errorbar(aes(year, ymin=ci_l, ymax=ci_u), data=obs, width=0, col=data_color) +
    geom_point(aes(year, prev, group=site, shape=type), data=ancsiteobs, col="darkgreen", alpha=0.5, size=0.8) +
    geom_line(aes(year, prev, group=site), data=ancsiteobs, col="darkgreen", alpha=0.5, size=0.3) +
    geom_line(aes(year, mean, color=model)) +
    geom_ribbon(aes(year, ymin=lower, ymax=upper, fill=model), alpha=0.3) +
    facet_wrap(~agegr) + ylim(0, max(df$upper)) +
    th +
    theme(legend.position="bottom",
          plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "in")) +
    ggtitle(paste0(out$country, " - ", out$eppregion), "Pregnant women prevalence")

  return(list(ggCore, ggAgegr3, ggPregPrev))
}
