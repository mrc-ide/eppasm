# HIV CLASS
# -----------------------------------------------------------------------------
hivEPP <- R6::R6Class("hivepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F, 
    public = list(
        year = 1,
        MODEL = "integer",
        data = "array",
        data_db = "array",
        grad = "array",
        grad_db = "array",
        f_death = "array",
        f_death_db = "array",
        initialize = function(fp, MODEL) {
            super$initialize(fp)
            MODEL <<- MODEL
            data <<- array(0, c(hDS, hAG, NG, PROJ_YEARS))
            f_death <<- grad <<- array(0, c(hDS, hAG, NG))
            if (MODEL==2) {
                data_db <<- data
                grad_db <<- grad
                f_death_db <<- f_death
            }
        })
)

hivFunc <- c(
aging = function(ag_prob) {
    data[,,,year]     <<- data[,,,year-1]
    nHup <- sweep(data[,-hAG,, year], 2:3, ag_prob[-hAG,], "*")
    data[,-hAG,,year] <<- data[,-hAG,,year] - nHup
    data[,  -1,,year] <<- data[,  -1,,year] + nHup
    if (MODEL==2) {
        data_db[,,,year]     <<- data_db[,,,year-1]
        nHup <- sweep(data_db[,-hAG,, year], 2:3, ag_prob[-hAG,], "*")
        data_db[,-hAG,,year] <<- data_db[,-hAG,,year] - nHup
        data_db[,  -1,,year] <<- data_db[,  -1,,year] + nHup
    }
},

add_entrants = function(artYesNo) {
    artNO = tail(artYesNo,2)
    if (MODEL==1)
        data[,1,,year] <<- data[,1,,year] + 
            sweep(p$paedsurv_cd4dist[,,year], 2, artNO, "*")
    if (MODEL==2) # add to virgin then debut
        data_db[,1,,year] <<- data_db[,1,,year] + 
            sweep(p$paedsurv_cd4dist[,,year], 2, artNO, "*")
},

sexual_debut = function() {
    debut_hiv <- sweep(data_db[,db_agr,, year], 2:3, p$db_pr, '*')
       data[,db_agr,,year] <<-    data[,db_agr,,year] + debut_hiv
    data_db[,db_agr,,year] <<- data_db[,db_agr,,year] - debut_hiv
},

deaths = function(survival_pr) {
    data[,,,year] <<- sweep(data[,,,year], 2:3, survival_pr, "*")
    if (MODEL==2)
        data_db[,,,year] <<- sweep(data_db[,,,year], 2:3, survival_pr, "*")
},

migration = function(migration_pr) {
    data[,,,year] <<- sweep(data[,,,year], 2:3, migration_pr, "*")
    if (MODEL==2)
        data_db[,,,year] <<- sweep(data_db[,,,year], 2:3, migration_pr, "*")
},

update_infection = function(new_infect) {
    grad[,,] <<- 0 # reset every time step
    grad <<- grad + sweep(p$cd4_initdist, 2:3, sumByAGs(new_infect, ag.idx), "*")
},

grad_progress = function(mortality_rate) { # HIV gradient progress
    if (p$eppmod == "directincid")
        grad[,,] <<- 0 # reset every time step
    # remove cd4 stage progression (untreated)
    nARTup <- p$cd4_prog * data[-hDS,,,year]
    grad[-hDS,,] <<- grad[-hDS,,] - nARTup
    grad[-1,,]   <<- grad[-1,,]   + nARTup # add 
    f_death <<- mortality_rate * data[,,,year]
    grad <<- grad - f_death # HIV mortality, untreated
    if (MODEL==2) {
        grad_db[,,] <<- 0 # reset every time, this's the 1st time grad_db is used
        nARTup <- p$cd4_prog * data_db[-hDS,,,year]
        grad_db[-hDS,,] <<- grad_db[-hDS,,] - nARTup
        grad_db[-1,,]   <<- grad_db[-1,,]   + nARTup 
        f_death_db <<- mortality_rate * data_db[,,,year]
        grad_db <<- grad_db - f_death_db
    }
},

eligible_for_art = function() {
    1 - (1 - rep(0:1, times=c(      p$artcd4elig_idx[year]-1,
                              hDS - p$artcd4elig_idx[year]+1))) * 
        (1 - rep(c(0, p$who34percelig),
                 c(2, hDS-2))) * 
        (1 - rep(p$specpop_percelig[year], hDS))
},

distribute_artinit = function(artinit, artpop) {
    debut_now     <- data_db[,,,year] + DT * grad_db
    all_hivpop    <- (data[,,,year] + DT * grad) + debut_now
    artinit       <- pmin(artinit, all_hivpop)
    pr.weight_db  <- debut_now / all_hivpop
    artinit_db    <- artinit * pr.weight_db
    artinit       <- artinit - artinit_db
    grad_db      <<- grad_db - artinit_db / DT
    artpop$grad_db_init(artinit_db)
    artinit
},

add_grad_to_pop = function() {
    data[,,,year] <<- data[,,,year] + DT * grad
    if (MODEL==2)
        data_db[,,,year] <<- data_db[,,,year] + DT * grad_db
},

adjust_pop = function(adj_prob) {
    data[,,,year] <<- sweep(data[,,,year], 2:3, adj_prob, "*")
    if (MODEL==2)
        data_db[,,,year] <<- sweep(data_db[,,,year], 2:3, adj_prob, "*")
},

# avoid default set method conflict
set_data = function(FUN="+", x, DS=T, AG=T, NG=T, YEAR=NULL) {
    FUN <- match.fun(FUN)
    if (is.null(YEAR))
        data[DS, AG, NG] <<- FUN(data[DS, AG, NG], x)
    else 
        data[DS, AG, NG, YEAR] <<- FUN(data[DS, AG, NG, YEAR], x)
}, 

get = function(YEAR=NULL, DS=T, AG=T, NG=T) {
    'get does not change anything, just return the data, use set to change'
    if (is.null(YEAR))
        data[DS, AG, NG, TRUE] # get all years
    else 
        data[DS, AG, NG, YEAR]
}, 

sweep_sex = function(FUN="*", x, year) {
    data[,,,year] <<- sweep(data[,,,year], 2:3, x, FUN)
})

setMembers(hivEPP, "public", names(hivFunc), hivFunc)