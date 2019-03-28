# EPP parameter super class
EPPClass <- setRefClass("eppfp", list(par="list", ss="list"))
EPPClass$methods(initialize = function(fp) {
    par <<- fp[-1]
    ss  <<- fp$ss
})

# POP CLASS
# -----------------------------------------------------------------------------
POPClass <- setRefClass("pop", list(data              = "array",
                                    MODEL             = "numeric", 
                                    birth_age         = "numeric", 
                                    birth_agrp        = "numeric", 
                                    prev15to49        = "vector", 
                                    incid15to49       = "vector", 
                                    entrantprev       = "vector", 
                                    pregprevlag       = "vector", 
                                    incrate15to49_ts  = "vector", 
                                    prev15to49_ts     = "vector", 
                                    birthslag         = "array", 
                                    infections        = "array", 
                                    hivdeaths         = "array", 
                                    natdeaths         = "array", 
                                    popadjust         = "array", 
                                    hivp_entrants_out = "array", 
                                    pop_db            = "array"), 
                                    contains          = "eppfp")

POPClass$methods(initialize = function(fp, MODEL=1) {
    'Initialize pop array'
    callSuper(fp)
    MODEL       <<- MODEL
    data        <<- array(0, c(ss$pAG, ss$NG, ss$pDS, ss$PROJ_YEARS))
    data[,,1,1] <<- par$basepop

    # Outputs
    prev15to49 <<- incid15to49 <<- entrantprev <<- pregprevlag <<- numeric(ss$PROJ_YEARS)
    infections <<- hivdeaths <<- natdeaths <<- popadjust <<- array(0, c(ss$pAG, ss$NG, ss$PROJ_YEARS))
    hivp_entrants_out    <<- array(0, c(ss$NG, ss$PROJ_YEARS))
    incrate15to49_ts     <<- rep(NA, length(par$rvec))
    prev15to49_ts        <<- rep(NA, length(par$rvec))
    birthslag            <<- par$birthslag

    if (MODEL==2) {
        pop_db       <<- array(0, c(ss$pDB, ss$NG, ss$pDS, ss$PROJ_YEARS))
        pop_db[,,,1] <<- sweep(data[1:ss$pDB,,,1], 1:2, 1 - par$db_pr, '*')
    }
})

POPClass$methods(
aging = function(year) {
    data[-c(1,ss$pAG),,,year] <<- data[-(ss$pAG-1:0),,,year-1]
    # open age group
    data[ss$pAG,,,year] <<- data[ss$pAG,,,year-1] + data[ss$pAG-1,,,year-1] 
}, 

set = function(FUN="+", x, AG=TRUE, NG=TRUE, DS=TRUE, YEAR=TRUE) {
    FUN <- match.fun(FUN)
    data[AG, NG, DS, YEAR] <<- FUN(data[AG, NG, DS, YEAR], x)
}, 

put = function(x, AG=TRUE, NG=TRUE, DS=TRUE, YEAR=TRUE) {
    'This can be done directly with `<-` using pop_object$data'
    data[AG, NG, DS, YEAR] <<- x
}, 

get = function(AG=TRUE, NG=TRUE, DS=TRUE, YEAR=TRUE) {
    'This can be done directly with `[` using pop_object$data'
    data[AG, NG, DS, YEAR]
}, 

sweep_sex = function(FUN="*", x, year) {
    data[,,,year] <<- sweep(data[,,,year], 1:2, x, FUN)
},

save_prev_n_inc = function(MODEL, year) {
    sus_pop <- data
    if (MODEL==2) 
      sus_pop[ss$db_aid,,,] <- sus_pop[ss$db_aid,,,] - pop_db
    prev15to49[year] <<- sum(data[ss$p.age15to49.idx,,ss$hivp.idx,year]) / 
                         sum(data[ss$p.age15to49.idx,,,year])
    incid15to49[year] <<- sum(incid15to49[year]) / 
                          sum(sus_pop[ss$p.age15to49.idx,,ss$hivn.idx,year-1])
},

update_fertile = function(MODEL, year, fp) {
    sus_pop <- data
    if (MODEL==2) 
      sus_pop[ss$db_aid,,,] <- sus_pop[ss$db_aid,,,] - pop_db
    birth_age <<- rowSums(sus_pop[ss$p.fert.idx,ss$f.idx,,year-1:0]) / 2 * par$asfr[,year]
    birth_agrp <<- sumByAG(birth_age, fp, TRUE)
    births     <- par$srb[,year] * sum(birth_agrp)
    if (year + ss$AGE_START <= ss$PROJ_YEARS)
      birthslag[, year + ss$AGE_START - 1] <<- births
}, 

cal_prev_pregant = function(MODEL, year, fp, hivpop, artpop) {
    sus_pop <- data
    if (MODEL==2) 
        sus_pop[ss$db_aid,,,] <- sus_pop[ss$db_aid,,,] - pop_db
    years <- year - 1:0
    meanWomen <- rowMeans(sus_pop[ss$p.fert.idx, ss$f.idx, ss$hivn.idx, years])
    hivn <- sumByAG(meanWomen, fp, TRUE)
    hivp <- rowMeans(hivpop$get(AG=ss$h.fert.idx, NG=ss$f.idx, YEAR=years),,2)
    art <- rowMeans(artpop$get(AG=ss$h.fert.idx, NG=ss$f.idx, YEAR=years),,3)
    pregprev <- sum(birth_agrp * (1 - hivn / (hivn + colSums(par$frr_cd4[,, year] * hivp) + colSums(par$frr_art[,,,year] * art,,2)))) / sum(birth_age)
    pregprevlag[year + ss$AGE_START - 1] <<- pregprev
})

# HIV CLASS
# -----------------------------------------------------------------------------
HIVClass <- setRefClass("hivpop", list(data="array"), contains="eppfp")
HIVClass$methods(initialize = function(fp, single_year = TRUE) {
    'Initialize hiv array'
    callSuper(fp)
    if (single_year)
        data <<- array(0, c(ss$hDS, ss$hAG, ss$NG))
    else
        data <<- array(0, c(ss$hDS, ss$hAG, ss$NG, ss$PROJ_YEARS))
})

# HIV gradient progress
HIVClass$methods(

progress = function(hiv_pop, mx) {
    # remove cd4 stage progression (untreated)
    data[-ss$hDS,,] <<- data[-ss$hDS,,] - par$cd4_prog * hiv_pop[-ss$hDS,,]
    # add cd4 stage progression (untreated)
    data[-1,,] <<- data[-1,,] + par$cd4_prog * hiv_pop[-ss$hDS,,]
    # HIV mortality, untreated
    data <<- data - mx * hiv_pop[,,]
}, 

# Age hiv 
aging = function(ag_prob, noart, year) {
    data[,,,year] <<- data[,,,year-1]
    nHup <- sweep(data[,-ss$hAG,, year], 2:3, ag_prob[-ss$hAG,], "*")
    data[,-ss$hAG,,year] <<- data[,-ss$hAG,,year] - nHup
    data[,-1,,year] <<- data[,-1,,year] + nHup
    data[, 1,,year] <<- data[,1,,year] + sweep(par$paedsurv_cd4dist[,,year], 2, noart, "*")
}, 

set = function(FUN="+", x, DS=TRUE, AG=TRUE, NG=TRUE, YEAR=NULL) {
    FUN <- match.fun(FUN)
    if (is.null(YEAR))
        data[DS, AG, NG] <<- FUN(data[DS, AG, NG], x)
    else 
        data[DS, AG, NG, YEAR] <<- FUN(data[DS, AG, NG, YEAR], x)
}, 

get = function(DS=TRUE, AG=TRUE, NG=TRUE, YEAR=NULL) {
    'get does not change anything, just return the data, use set to change'
    if (is.null(YEAR))
        data[DS, AG, NG]
    else 
        data[DS, AG, NG, YEAR]
}, 

sweep_sex = function(FUN="*", x, year) {
    data[,,,year] <<- sweep(data[,,,year], 2:3, x, FUN)
})

# ART GRAD
# -----------------------------------------------------------------------------
ARTClass <- setRefClass("artpop", list(data="array"), contains="eppfp")

ARTClass$methods(initialize = function(fp, single_year=TRUE) {
    'Initialize hiv array'
    callSuper(fp)
    if (single_year)
        data <<- array(0, c(ss$hTS, ss$hDS, ss$hAG, ss$NG))
    else
        data <<- array(0, c(ss$hTS, ss$hDS, ss$hAG, ss$NG, ss$PROJ_YEARS))
})

ARTClass$methods(
## progression and mortality (HARD CODED 6 months duration)
progress = function(art_pop, year) {
    # remove ART duration progression 
    data[1:2,,,] <<- data[1:2,,,] - 2.0 * art_pop[1:2,,,]
    # add ART duration progression
    data[2:3,,,] <<- data[2:3,,,] + 2.0 * art_pop[1:2,,,] 
    # ART mortality
    data <<- data - par$art_mort * par$artmx_timerr[,year] * art_pop
},

aging = function(ag_prob, isart, year) {
    data[,,,,year] <<- data[,,,,year-1]
    nARTup <- sweep(data[,,-ss$hAG,,year], 3:4, ag_prob[-ss$hAG,], "*")
    data[,,-ss$hAG,,year] <<- data[,,-ss$hAG,,year] - nARTup
    data[,,-1,,year] <<- data[,,-1,,year] + nARTup
    data[,,1,,year] <<- data[,,1,,year] + 
                        sweep(par$paedsurv_artcd4dist[,,,year], 3, isart, "*")
},

set = function(FUN="+", x, TS=TRUE, DS=TRUE, AG=TRUE, NG=TRUE, YEAR=NULL) {
    FUN <- match.fun(FUN)
    if (is.null(YEAR))
        data[TS, DS, AG, NG] <<- FUN(data[TS, DS, AG, NG], x)
    else
        data[TS, DS, AG, NG, YEAR] <<- FUN(data[TS, DS, AG, NG, YEAR], x)
}, 

get = function(TS=TRUE, DS=TRUE, AG=TRUE, NG=TRUE, YEAR=NULL) {
    'get does not change anything, just return the data, use set to change'
    if (is.null(YEAR))
        data[TS, DS, AG, NG]
    else
        data[TS, DS, AG, NG, YEAR]
},

sweep_sex = function(FUN="*", x, year) {
    data[,,,,year] <<- sweep(data[,,,,year], 3:4, x, FUN)
})