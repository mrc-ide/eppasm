#' Annualized number of new infections
#'
calc_infections_eppspectrum <- function(fp, pop, hivpop, i, ii, r_ts){

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

  art.ii <- sum(hivpop[-1,,h.age15to49.idx,,i])
  if(sum(hivpop[,,h.age15to49.idx[1],,i]) > 0)
    art.ii <- art.ii - sum(pop[p.age15to49.idx[1],,hivp.idx,i] * colSums(hivpop[-1,,h.age15to49.idx[1],,i],,2) / colSums(hivpop[,,h.age15to49.idx[1],,i],,2)) * (1-DT*(ii-1))
  if(sum(hivpop[,,tail(h.age15to49.idx, 1)+1,,i]) > 0)
    art.ii <- art.ii + sum(pop[tail(p.age15to49.idx,1)+1,,hivp.idx,i] * colSums(hivpop[-1,,tail(h.age15to49.idx, 1)+1,,i],,2) / colSums(hivpop[,,tail(h.age15to49.idx, 1)+1,,i],,2)) * (1-DT*(ii-1))
  
  transm_prev <- (hivp.ii - art.ii + fp$relinfectART*art.ii) / (hivn.ii+hivp.ii)

  incrate15to49.ts <- r_ts * transm_prev + fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
  sexinc15to49.ts <- incrate15to49.ts*c(1, fp$incrr_sex[i])*sum(pop[p.age15to49.idx,,hivn.idx,i])/(sum(pop[p.age15to49.idx,m.idx,hivn.idx,i]) + fp$incrr_sex[i]*sum(pop[p.age15to49.idx, f.idx,hivn.idx,i]))
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
  infections.ts <- agesex.inc * pop[,,hivn.idx,i]

  attr(infections.ts, "incrate15to49.ts") <- incrate15to49.ts
  attr(infections.ts, "prevcurr") <- hivp.ii / (hivn.ii+hivp.ii)

  return(infections.ts)
}

calc_infections_simpletransm <- function(fp, pop, hivpop, i, ii, r_ts){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  ts <- (i-2)/DT + ii

  ## Calculate prevalence of unsuppressed viral load among sexually active population
  hivn.ii <- colSums(pop[p.age15to49.idx,,hivn.idx,i])
  hivn.ii <- hivn.ii - pop[p.age15to49.idx[1],,hivn.idx,i]*(1-DT*(ii-1))
  hivn.ii <- hivn.ii + pop[tail(p.age15to49.idx,1)+1,,hivn.idx,i]*(1-DT*(ii-1))

  ## Calculate proportion in each HIV age group who are in 15 to 49 population, accounting for partial year time step
  ha1 <- h.age15to49.idx[1]  
  haM <- h.age15to49.idx[length(h.age15to49.idx)]+1  # age group one above 15 to 49
  prop_include <- rbind(ifelse(pop[agfirst.idx[ha1],,hivp.idx,i] > 0,
                               1 - pop[agfirst.idx[ha1],,hivp.idx,i] / colSums(pop[agfirst.idx[ha1]+1:h.ag.span[ha1]-1,,hivp.idx,i]) * (1-DT*(ii-1)),
                               c(1.0, 1.0)),
                        matrix(1, length(h.age15to49.idx)-1, NG),
                        ifelse(pop[agfirst.idx[haM],,hivp.idx,i] > 0,
                               pop[agfirst.idx[haM],,hivp.idx,i] / colSums(pop[agfirst.idx[haM]+1:h.ag.span[haM]-1,,hivp.idx,i]) * (1-DT*(ii-1)),
                               c(0, 0)))
                        

  hivp_noart.ii <- colSums(colSums(sweep(hivpop[1,,c(h.age15to49.idx, haM),,i], 1, fp$relsexact_cd4cat, "*")) * prop_include)
  art.ii <- colSums(colSums(hivpop[-1,,c(h.age15to49.idx, haM),,i],,2) * prop_include)
  
  ## Prevalence of unsuppressed viral load among sexually active population
  hivtransm_prev <- (hivp_noart.ii + fp$relinfectART * art.ii) / (hivn.ii+hivp_noart.ii+art.ii)
  
  ## r_sex[1:2] is the transmission rate by (Men, Women)
  r_sex <- c(sqrt(fp$mf_transm_rr), 1/sqrt(fp$mf_transm_rr)) * r_ts

  sexinc15to49.ts <- (r_sex * hivtransm_prev)[2:1] + fp$mf_transm_rr^c(-0.25, 0.25) * fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")
  infections.ts <- agesex.inc * pop[,,hivn.idx,i]

  attr(infections.ts, "incrate15to49.ts") <- sum(infections.ts[p.age15to49.idx,]) / sum(hivn.ii)
  attr(infections.ts, "prevcurr") <- sum(hivp_noart.ii+art.ii) / sum(hivn.ii+hivp_noart.ii+art.ii)

  return(infections.ts)
}
