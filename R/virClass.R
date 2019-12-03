# Virgin CLASS
# -----------------------------------------------------------------------------
virginEPP <- R6::R6Class("virginepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F, 
    public = list(
        year = 1,
        MODEL = "integer",
        data = "array",
        initialize = function(fp, MODEL) {                               
            super$initialize(fp)
            MODEL <<- MODEL
            data  <<- array(0, c(pDB, NG, pDS, PROJ_YEARS)) # debut ages only
            data[,,hivn.idx,1] <<- p$basepop[db_aid,] * (1 - p$db_pr)
        })
)

virginFunc <- c(
aging = function() {
    data[-1,,,year] <<- data[-pDB,,,year-1]
},

sexual_debut = function() {
    debut_now      <- sweep(data[,,,year], 1:2, p$db_pr, "*")
    data[,,,year] <<- data[,,,year] - debut_now
},

deaths = function() {
    death <- sweep(data[,,,year], 1:2, 1 - p$Sx[db_aid,,year], "*")
    data[,,,year] <<- data[,,,year] - death
},

remove_hiv_death = function(p_hiv_death_by_age) {
    data[,,hivp.idx,year] <<- data[,,hivp.idx,year] * (1 - p_hiv_death_by_age[db_aid,])
},

migration = function(migration_pr) {
    data[,,,year] <<- sweep(data[,,,year], 1:2, migration_pr[db_aid,], "*")
},

adjust_pop = function(adj_prob) {
    data[,,,year] <<- sweep(data[,,,year], 1:2, adj_prob[db_aid,,year], "*")
}, 

n_all = function(dt=0, only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    sum(data[age_range,,,year]) -
        sum(data[age_range[1],,,year]) * dt +
            sum(data[tail(age_range)+1,,,year]) * dt
},

n_HIV = function(dt=0, only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    sum(data[age_range,, hivp.idx, year]) -
        sum(data[age_range[1],, hivp.idx, year]) * dt +
            sum(data[tail(age_range)+1,, hivp.idx, year]) * dt
},

n_NEG = function(only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    data[age_range,, hivn.idx, year]
},

n_IDU = function(only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    data[age_range,,,year] * p$idu_virgin_pro[age_range,,]
}
)

setMembers(virginEPP, "public", names(virginFunc), virginFunc)