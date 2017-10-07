#' Annualized number of new infections
#'
calc_infections_eppspectrum <- function(fp, pop, hivpop, i, ii, prevlast, rveclast){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  ## HIV population size at ts
  ts <- (i-2)/DT + ii

  hivn.ii <- sum(pop[p.age15to49.idx,,hivn.idx,i])
  hivn.ii <- hivn.ii - sum(pop[p.age15to49.idx[1],,hivn.idx,i])*(1-DT*(ii-1))
  hivn.ii <- hivn.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivn.idx,i])*(1-DT*(ii-1))

  hivp.ii <- sum(pop[p.age15to49.idx,,hivp.idx,i])
  hivp.ii <- hivp.ii - sum(pop[p.age15to49.idx[1],,hivp.idx,i])*(1-DT*(ii-1))
  hivp.ii <- hivp.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivp.idx,i])*(1-DT*(ii-1))

  ## there is an approximation here since this is the 15-49 pop (doesn't account for the slight offset in age group)
  propart.ii <- ifelse(hivp.ii > 0, sum(hivpop[-1,,h.age15to49.idx,,i])/sum(hivpop[,,h.age15to49.idx,,i]), 0)

  ## incidence

  ## calculate r(t)
  prevcurr <- hivp.ii / (hivn.ii+hivp.ii)
  if(fp$eppmod %in% c("rtrend", "rtrend_rw"))
    r_ts <- calc_rtrend_rt(fp$proj.steps[ts], fp, rveclast, prevlast, prevcurr)
  else
    r_ts <- fp$rvec[ts]

  transm_prev <- hivp.ii * (1 - (1-fp$relinfectART)*propart.ii) / (hivn.ii+hivp.ii)

  incrate15to49.ts <- r_ts * hivp.ii * (1 - (1-fp$relinfectART)*propart.ii) / (hivn.ii+hivp.ii) + fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
  sexinc15to49.ts <- incrate15to49.ts*c(1, fp$incrr_sex[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$incrr_sex[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
  infections.ts <- agesex.inc * pop[,,hivn.idx,i]

  attr(infections.ts, "rvec") <- r_ts
  attr(infections.ts, "incrate15to49.ts") <- incrate15to49.ts
  attr(infections.ts, "prevcurr") <- incrate15to49.ts

  return(infections.ts)
}


calc_rtrend_rt <- function(t, fp, rveclast, prevlast, prevcurr){
  if(t > fp$tsEpidemicStart){
    par <- fp$rtrend
    gamma.t <- if(t < par$tStabilize) 0 else (prevcurr-prevlast)*(t - par$tStabilize) / (fp$ss$DT*prevlast)
    logr.diff <- par$beta[2]*(par$beta[1] - rveclast) + par$beta[3]*prevlast + par$beta[4]*gamma.t
      return(exp(log(rveclast) + logr.diff))
    } else
      return(fp$rtrend$r0)
}

