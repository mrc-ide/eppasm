#' Annualized number of new infections
#'
calc_infections_eppspectrum <- function(fp, mx, pop, hivpop, artpop,
                                        i, ii, r_ts, isMixing, FOI){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  if (!isMixing) {
  ## HIV population size at ts
    ts <- (i-2)/DT + ii
    wDT <- 1-DT*(ii-1)

    nHts <- function(y, p)
      sum(pop[p,,y,i]) - sum(pop[p[1],,y,i]) * wDT + sum(pop[tail(p,1)+1,,y,i]) * wDT

    nARTa <- function(p, h) {
      num <- colSums(artpop[,,h,,i],,2)
      den <- num + colSums(hivpop[,h,,i],,1)
      sum(pop[p,,hivp.idx,i] * num / den)
    }

    nARTts <- function(p, h) {
      art.ii <- sum(artpop[,,h,,i])
      if (sum(hivpop[,h[1],,i]) + sum(artpop[,,h[1],,i])  > 0)
        art.ii <- art.ii - nARTa(p[1], h[1]) *  wDT
      if (sum(hivpop[,tail(h, 1)+1,,i]) + sum(artpop[,,tail(h, 1)+1,,i]) > 0)
        art.ii <- art.ii + nARTa(tail(p,1)+1, tail(h, 1)+1) *  wDT
      return(art.ii)
    }

    hivn.ii <- nHts(hivn.idx, p.age15to49.idx)
    hivp.ii <- nHts(hivp.idx, p.age15to49.idx)
    art.ii  <- nARTts(p.age15to49.idx, h.age15to49.idx)

    # TODO: not fixing for number of risk groups
    transm_prev <- (hivp.ii - art.ii + fp$relinfectART*art.ii) / (hivn.ii+hivp.ii)
    incrate15to49.ts <- r_ts * transm_prev + fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
    tmp1_ <- function(p, y) sum(pop[p,,y,i]) / (sum(pop[p,m.idx,y,i]) + fp$incrr_sex[i]*sum(pop[p, f.idx,y,i]))
    sexinc15to49.ts <- incrate15to49.ts * c(1, fp$incrr_sex[i]) * tmp1_(p.age15to49.idx, hivn.idx)
    tmp2_ <- function(p, y) colSums(sweep(pop[p,,y,,i], 1:2, fp$incrr_age[p,,i], '*')) / colSums(pop[p,,y,,i])
    agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts / rowSums(tmp2_(p.age15to49.idx, hivn.idx)), "*")
    infections.ts <- sweep(pop[,,hivn.idx,,i], 1:2, agesex.inc, '*')
  } else {
    FOI[,,i] <- FOIs(pop, i, mx, fp)
    infections.ts <- pop[,,hivn.idx,i] * FOI[,,i]
    if(any(infections.ts < 0) | any(infections.ts > 1000)) browser()
    incrate15to49.ts <- NA # TODO: ?
  }

  attr(infections.ts, "incrate15to49.ts") <- incrate15to49.ts
  attr(infections.ts, "prevcurr") <- NA
  return(infections.ts)
}

## Beers coefficients to distribute infections from 5-year age groups to single-year of age
create_beers <- function(n5yr){

  ## Beer's coefficients for disaggregating 5 year age groups into
  ## single-year age groups (from John Stover)
  Afirst <- rbind(c(0.3333, -0.1636, -0.0210,  0.0796, -0.0283),
                  c(0.2595, -0.0780,  0.0130,  0.0100, -0.0045),
                  c(0.1924,  0.0064,  0.0184, -0.0256,  0.0084),
                  c(0.1329,  0.0844,  0.0054, -0.0356,  0.0129),
                  c(0.0819,  0.1508, -0.0158, -0.0284,  0.0115))
  Asecond <- rbind(c( 0.0404,  0.2000, -0.0344, -0.0128,  0.0068),
                   c( 0.0093,  0.2268, -0.0402,  0.0028,  0.0013),
                   c(-0.0108,  0.2272, -0.0248,  0.0112, -0.0028),
                   c(-0.0198,  0.1992,  0.0172,  0.0072, -0.0038),
                   c(-0.0191,  0.1468,  0.0822, -0.0084, -0.0015))
  Amid <- rbind(c(-0.0117,  0.0804,  0.1570, -0.0284,  0.0027),
                c(-0.0020,  0.0160,  0.2200, -0.0400,  0.0060),
                c( 0.0050, -0.0280,  0.2460, -0.0280,  0.0050),
                c( 0.0060, -0.0400,  0.2200,  0.0160, -0.0020),
                c( 0.0027, -0.0284,  0.1570,  0.0804, -0.0117))
  Apenult <- rbind(c(-0.0015, -0.0084,  0.0822,  0.1468, -0.0191),
                   c(-0.0038,  0.0072,  0.0172,  0.1992, -0.0198),
                   c(-0.0028,  0.0112, -0.0248,  0.2272, -0.0108),
                   c( 0.0013,  0.0028, -0.0402,  0.2268,  0.0093),
                   c( 0.0068, -0.0128, -0.0344,  0.2000,  0.0404))
  Aultim <- rbind(c( 0.0115, -0.0284, -0.0158,  0.1508,  0.0819),
                  c( 0.0129, -0.0356,  0.0054,  0.0844,  0.1329),
                  c( 0.0084, -0.0256,  0.0184,  0.0064,  0.1924),
                  c(-0.0045,  0.0100,  0.0130, -0.0780,  0.2595),
                  c(-0.0283,  0.0796, -0.0210, -0.1636,  0.3333))

  A <- do.call(rbind,
               c(list(cbind(Afirst, matrix(0, 5, n5yr-5)),
                      cbind(Asecond, matrix(0, 5, n5yr-5))),
                 lapply(0:(n5yr-6), function(i) cbind(matrix(0, 5, i), Amid, matrix(0, 5, (n5yr-5)-i))),
                 list(cbind(matrix(0, 5, n5yr-6), Apenult, matrix(0, 5, 1)),
                      cbind(matrix(0, 5, n5yr-6), Aultim, matrix(0, 5, 1)),
                      c(rep(0, n5yr-1), 1))))
  return(round(A, 4))
}
  
## Beers coefficient matrix
beers_Amat <- create_beers(17)[16:81, 4:17]

