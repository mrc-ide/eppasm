
#' @import data.table
#' 
tidy_output <- function(fit, modlab, country=NA, eppregion=NA, ancsite = TRUE){

  idvars <- data.frame(country = country,
                       eppregion = eppregion,
                       modlab = modlab)
  print(paste(country, eppregion))

  ss <- fit$fp$ss

  
  ## agegr3 predictions
  agegr3 <- merge(expand.grid(year = 1999:2018, sex=c("male", "female"),
                              age = c(15, 25, 35)),
                  data.frame(age = c(15, 25, 35), agspan = c(10, 10, 15)))
  agegr3$agegr3 <- with(agegr3, paste0(age, "-", age+agspan-1))
  agegr3$sidx <- match(agegr3$sex, c("both", "male", "female")) - 1L
  agegr3$aidx <- agegr3$age - ss$AGE_START + 1
  agegr3$yidx <- agegr3$year - ss$proj_start + 1


  ## simulate model projections
  param_list <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))
  
  if(fit$fp$eppmod == "rspline")
    fit <- rw_projection(fit)
  
  fp_list <- lapply(param_list, function(par) update(fit$fp, list=par))
  mod_list <- lapply(fp_list, simmod)



  ## Assemble results

  year <- fit$fp$ss$proj_start + 1:fit$fp$ss$PROJ_YEARS - 1L


  ## Core results
  rvec <- sapply(mod_list, attr, "rvec_ts")
  rvecidx <- c(1, seq_len(ss$PROJ_YEARS-1) * ss$hiv_steps_per_year)
  rvec <- rvec[rvecidx, ]

  prev <- sapply(mod_list, prev)
  incid <- mapply(incid, mod = mod_list, fp = fp_list)
  popsize <- sapply(mod_list, colSums, dims=3)
  plhiv <- sapply(lapply(mod_list, function(x) x[,,2,]), colSums, dims = 2)
  infections <- sapply(lapply(mod_list, attr, "infections"), colSums, dims=2)
  
  aidsdeaths <- sapply(lapply(mod_list, attr, "hivdeaths"), colSums, dims=2)
  artcov <- sapply(mod_list, artcov15plus)
  ancartcov <- mapply(agepregartcov, mod_list, fp_list, MoreArgs=list(aidx=1, yidx=1:fit$fp$ss$PROJ_YEARS, agspan=35, expand=TRUE))

  artcov[is.na(artcov)] <- 0
  ancartcov[is.na(ancartcov)] <- 0

  year <- fit$fp$ss$proj_start + 1:fit$fp$ss$PROJ_YEARS - 1L
  
  core <- list(prev = prev,
               incid = incid,
               aidsdeaths15pl = aidsdeaths,
               artcov15pl = artcov,
               pregartcov = ancartcov,
               popsize15pl = popsize,
               plhiv15pl = plhiv,
               infections15pl = infections,
               "r(t)" = rvec,
               "log r(t)" = log(rvec))
  core <- lapply(core, estci2)
  core <- Map(data.frame,
              list(idvars),
              year = list(year),
              indicator = names(core),
              core)
  core <- do.call(rbind, core)

  ## ADD INCIDENCE SEX RATIO
               
  ## core <- rbind(core,
  ##               data.frame(country = country,
  ##                          eppregion = eppregion,
  ##                          year = dimnames(fit$agegr3incid)[[3]],
  ##                          model = modlab,
  ##                          indicator = "incidence sex ratio",
  ##                          estci2(fit$agegr3incid["15-49", "female",,] /
  ##                                 fit$agegr3incid["15-49", "male",,])))


  agegr3prev <- mapply(ageprev, mod=mod_list,
                       MoreArgs=list(aidx=agegr3$aidx,
                                     sidx=agegr3$sidx,
                                     yidx=agegr3$yidx,
                                     agspan=agegr3$agspan))
  
  agegr3prev <- data.frame(country = country,
                           eppregion = eppregion,
                           agegr3,
                           model = modlab,
                           estci2(agegr3prev))

  ## Prevalence among pregnant women
  pregprev_pred <- rbind(data.frame(year = 1980:2020, age = 15, agspan = 35),
                         data.frame(year = 1980:2020, age = 15, agspan = 10),
                         data.frame(year = 1980:2020, age = 25, agspan = 25),
                         expand.grid(year = 1980:2020, age = 3:9*5, agspan = 5))
  pregprev_pred$agegr <- with(pregprev_pred, paste0(age, "-", age+agspan-1))
  pregprev_pred$agegr <- factor(pregprev_pred$agegr, c("15-49", "15-24", "25-49", paste0(3:9*5, "-", 3:9*5+4)))
  pregprev_pred$aidx <- pregprev_pred$age - ss$AGE_START + 1
  pregprev_pred$yidx <- pregprev_pred$year - ss$proj_start+ 1

  pregprev <- mapply(agepregprev, mod=mod_list, fp=fp_list,
                     MoreArgs = list(aidx=pregprev_pred$aidx,
                                     yidx=pregprev_pred$yidx,
                                     agspan=pregprev_pred$agspan))
  pregprev <- data.frame(idvars,
                         pregprev_pred,
                         estci2(pregprev))

  
  ## agegr3incid <- reshape2::melt(estci2(fit$agegr3incid),
  ##                       varnames = c("agegr", "sex", "year", "outcome"))
  ## ageincid <- data.frame(country = country,
  ##                        eppregion = eppregion,
  ##                        indicator = "incidence",
  ##                        model=modlab,
  ##                        reshape2::dcast(ageincid, ... ~ outcome, value.var="value"))

  ## ageinfections <- reshape2::melt(estci2(fit$agegr3incid),
  ##                                 varnames = c("agegr", "sex", "year", "outcome"))
  ## ageinfections <- data.frame(country = country,
  ##                             eppregion = eppregion,
  ##                             indicator = "infections",
  ##                             model=modlab,
  ##                             reshape2::dcast(ageinfections, ... ~ outcome, value.var="value"))

  ## relincid <- reshape2::melt(fit$relincid,
  ##                            varnames = c("agegr", "sex", "year", "outcome"))
  ## relincid <- data.frame(country = country,
  ##                        eppregion = eppregion,
  ##                        indicator = "relincid",
  ##                        model=modlab,
  ##                        reshape2::dcast(relincid, ... ~ outcome, value.var="value"))

  ## Age specific prevalence predictions in survey years
  vnames <- c("year", "sex", "agegr", "n", "prev", "se", "ci_l", "ci_u", "aidx", "sidx", "yidx", "agspan")
  ageprevdat <- fit$likdat$hhs.dat[intersect(names(fit$likdat$hhs.dat), vnames)]

  ageprevpred <- mapply(ageprev,
                       mod = mod_list,
                       MoreArgs = list(aidx=ageprevdat$aidx,
                                       sidx=ageprevdat$sidx,
                                       yidx=ageprevdat$yidx,
                                       agspan=ageprevdat$agspan))
  ageprevdat <- data.frame(country = country,
                           eppregion = eppregion,
                           model = modlab,
                           ageprevdat,
                           estci2(ageprevpred))

  ## Site-level ANC outputs
  if(ancsite && nrow(fit$likdat$ancsite.dat$df)) {

    b_site <- Map(sample_b_site, mod_list, fp_list,
                  list(fit$likdat$ancsite.dat), resid = FALSE)

    b_site_sigma <- sapply(b_site, anclik::sample.sigma2)

    b_site_df <- estci2(do.call(cbind, b_site))
    ancsite_b <- data.frame(idvars, site = rownames(b_site_df), b_site_df)
    
    newdata <- expand.grid(site = unique(fit$likdat$ancsite.dat$df$site),
                           year = 1985:2020,
                           type = "ancss",
                           age = 15,
                         agspan = 35,
                         n = 300)
    new_df <- ancsite_pred_df(newdata, fit$fp)
    
    ancsite_pred <- mapply(sample_ancsite_pred, mod_list, fp_list,
                           b_site = b_site,
                           MoreArgs = list(newdata = new_df))
    ancsite_pred <- data.frame(newdata, estci2(ancsite_pred))
    
    ancsite_pred <- merge(ancsite_pred,
                          fit$likdat$ancsite.dat$df[c("site", "year", "type", "age", "agspan", "n", "prev", "pstar", "W", "v")],
                          by = c("site", "year", "type", "age", "agspan"),
                          suffixes = c("_sim", "_obs"), all.x=TRUE)
    
    ancsite_pred <- data.frame(idvars, ancsite_pred)
  } else {
    ancsite_pred <- NULL
    ancsite_b <- NULL
  }
   
  out <- list(core = core,
              ageprevdat=ageprevdat,
              pregprev = pregprev,
              agegr3prev = agegr3prev,
              ancsite_pred = ancsite_pred,
              ancsite_b = ancsite_b)
              ## ageincid = ageincid,
              ## ageinfections = ageinfections,
              ## relincid = relincid)
}


combine_tidy <- function(...){
  Map(rbind, ...)
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



get_pointwise_ll <- function(fit, newdata = fit$likdat){

  ## simulate model projections
  param_list <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))
  
  if(fit$fp$eppmod == "rspline")
    fit <- rw_projection(fit)
  if(fit$fp$eppmod == "rhybrid")
    fit <- extend_projection(fit, proj_years = fit$fp$ss$PROJ_YEARS)
  
  fp_list <- lapply(param_list, function(par) update(fit$fp, list=par))
  mod_list <- lapply(fp_list, simmod)


  ## Site-level ANC data
  if(nrow(fit$likdat$ancsite.dat$df)) {
    b_site <- Map(sample_b_site, mod_list, fp_list,
                  list(fit$likdat$ancsite.dat), resid = FALSE)
    
    ancsite_ll <- mapply(ll_ancsite_conditional, mod_list, fp_list,
                         b_site = b_site,
                         MoreArgs = list(newdata = newdata$ancsite.dat))
  } else {
    ancsite_ll <- NULL
  }

  ## ANC-RT census data
  if(exists("ancrtcens.dat", newdata))
    ancrtcens_ll <- mapply(ll_ancrtcens, mod = mod_list, fp = fp_list,
                           MoreArgs = list(dat = newdata$ancrtcens.dat, pointwise = TRUE))
  else
    ancrtcens_ll <- NULL

  ## Household survey data
  if(exists("hhs.dat", where=newdata)){
    if(exists("ageprev", fit$fp) && fit$fp$ageprev=="binom")
      hhs_ll <- mapply(ll_hhsage_binom, mod_list,
                       MoreArgs = list(dat = newdata$hhs.dat, pointwise = TRUE))
    else ## use probit likelihood
      hhs_ll <- mapply(ll_hhsage, mod_list,
                       MoreArgs = list(dat = newdata$hhs.dat, pointwise = TRUE))
  } else {
    hhs_ll <- NULL
  }

  list(ancsite = ancsite_ll,
       ancrtcens = ancrtcens_ll,
       hhs = hhs_ll)
}
