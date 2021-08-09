#' Aggregate a list of fits 
#' 
#' @details
#'
#' This makes some assumptions:
#' 
#' * The number of resamples is the same for all fits in list
#' * The state space dimension, including number of projection years, is the same
#' 
tidy_aggr <- function(fitlist, modlab = NA, geolab=NA){

  print(paste(geolab, modlab))

  fp_seq <- lapply(fitlist, "[[", "fp")
  resample_seq <- lapply(fitlist, "[[", "resample")
  
  ## simulate model projections
  aggr_one <- function(i) {

    param_seq <- Map(fnCreateParam,
                     lapply(resample_seq, function(x) x[i,]),
                     fp_seq)
    param_seq <- Map(update, fp_seq, list = param_seq)

    mod_seq <- lapply(param_seq, simmod)

    modaggr <- Reduce("+", mod_seq)
    attr(modaggr, "hivpop") <- Reduce("+", lapply(mod_seq, attr, "hivpop"))
    attr(modaggr, "artpop") <- Reduce("+", lapply(mod_seq, attr, "artpop"))
    attr(modaggr, "infections") <- Reduce("+", lapply(mod_seq, attr, "infections"))
    attr(modaggr, "hivdeaths") <- Reduce("+", lapply(mod_seq, attr, "hivdeaths"))
    attr(modaggr, "natdeaths") <- Reduce("+", lapply(mod_seq, attr, "natdeaths"))

    attr(modaggr, "prev15to49") <- calc_prev15to49(modaggr, param_seq[[1]])
    attr(modaggr, "incid15to49") <- calc_incid15to49(modaggr, param_seq[[1]])

    class(modaggr) <- c("specaggr", "spec")
    modaggr
  }

  mod_list <- lapply(seq_len(nrow(resample_seq[[1]])), aggr_one)
                     
  
  ## ID variables  
  idvars <- data.frame(geolab = geolab,
                       modlab = modlab)

  ss <- fitlist[[1]]$fp$ss
  year <- ss$proj_start + 1:ss$PROJ_YEARS - 1L
  
  ## agegr3 predictions
  agegr3 <- merge(expand.grid(year = 1999:2018, sex=c("male", "female"),
                              age = c(15, 25, 35)),
                  data.frame(age = c(15, 25, 35), agspan = c(10, 10, 15)))
  agegr3$agegr3 <- with(agegr3, paste0(age, "-", age+agspan-1))
  agegr3$sidx <- match(agegr3$sex, c("both", "male", "female")) - 1L
  agegr3$aidx <- agegr3$age - ss$AGE_START + 1
  agegr3$yidx <- agegr3$year - ss$proj_start + 1


  ## Core results

  prev <- sapply(mod_list, prev)
  incid <- mapply(incid, mod = mod_list, fp = list(NULL))
  
  aidsdeaths <- sapply(lapply(mod_list, attr, "hivdeaths"), colSums, dims=2)
  artcov <- sapply(mod_list, artcov15plus)
  artcov[is.na(artcov)] <- 0
  
  popsize <- sapply(mod_list, colSums, dims=3)
  plhiv <- sapply(lapply(mod_list, function(x) x[,,2,]), colSums, dims = 2)
  infections <- sapply(lapply(mod_list, attr, "infections"), colSums, dims=2)

  artcov15to49 <- sapply(mod_list, artcov15to49)
  artcov15to49[is.na(artcov15to49)] <- 0

  ## ASSUMES SAME relinfectART FOR ALL FITS
  rvec <- (incid + rbind(0, diff(incid))) / (prev * (1 - artcov15to49 * (1 - fp_seq[[1]]$relinfectART)))

  core <- list(prev = prev,
               incid = incid,
               aidsdeaths15pl = aidsdeaths,
               artcov15pl = artcov,
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
  
  agegr3prev <- data.frame(idvars,
                           agegr3,
                           estci2(agegr3prev))

  ## ## Prevalence among pregnant women
  ## pregprev_pred <- rbind(data.frame(year = 1980:2020, age = 15, agspan = 35),
  ##                        data.frame(year = 1980:2020, age = 15, agspan = 10),
  ##                        data.frame(year = 1980:2020, age = 25, agspan = 25),
  ##                        expand.grid(year = 1980:2020, age = 3:9*5, agspan = 5))
  ## pregprev_pred$agegr <- with(pregprev_pred, paste0(age, "-", age+agspan-1))
  ## pregprev_pred$agegr <- factor(pregprev_pred$agegr, c("15-49", "15-24", "25-49", paste0(3:9*5, "-", 3:9*5+4)))
  ## pregprev_pred$aidx <- pregprev_pred$age - ss$AGE_START + 1
  ## pregprev_pred$yidx <- pregprev_pred$year - ss$proj_start+ 1

  ## pregprev <- mapply(agepregprev, mod=mod_list, fp=fp_list,
  ##                    MoreArgs = list(aidx=pregprev_pred$aidx,
  ##                                    yidx=pregprev_pred$yidx,
  ##                                    agspan=pregprev_pred$agspan))
  ## pregprev <- data.frame(idvars,
  ##                        pregprev_pred,
  ##                        estci2(pregprev))

  
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
  ## vnames <- c("year", "sex", "agegr", "n", "prev", "se", "ci_l", "ci_u", "aidx", "sidx", "yidx", "agspan")
  ## ageprevdat <- fit$likdat$hhs.dat[intersect(names(fit$likdat$hhs.dat), vnames)]

  ## ageprevpred <- mapply(ageprev,
  ##                      mod = mod_list,
  ##                      MoreArgs = list(aidx=ageprevdat$aidx,
  ##                                      sidx=ageprevdat$sidx,
  ##                                      yidx=ageprevdat$yidx,
  ##                                      agspan=ageprevdat$agspan))
  ## ageprevdat <- data.frame(country = country,
  ##                          eppregion = eppregion,
  ##                          model = modlab,
  ##                          ageprevdat,
  ##                          estci2(ageprevpred))

  out <- list(core = core,
              ## ageprevdat=ageprevdat,
              ## pregprev = pregprev,
              agegr3prev = agegr3prev)

out
}
