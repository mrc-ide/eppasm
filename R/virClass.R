# Virgin CLASS
# -----------------------------------------------------------------------------
virginFunc <- c(
aging = function() {
    data[-1,,,year] <<- data[-pDB,,,year-1]
},

sexual_debut = function() {
    debut_now      <- sweep(data[,,,year], 1:2, p$db_rate[,,year], "*")
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
    sum(data[age_range,,,year]) - sum(data[age_range[1],,,year]) * dt
},

n_HIV = function(dt=0) {
    sum(data[db_aid,, hivp.idx, year]) -
        sum(data[db_aid[1],, hivp.idx, year]) * dt
},

n_NEG = function(dt=0) {
    data[db_aid,, hivn.idx, year] - 
        sum(data[db_aid[1],, hivp.idx, year]) * dt
},

n_IDU = function() {
    data[db_aid,,,year] * p$idu_virgin_pro[db_aid,,]
}
)

setMembers(virginEPP, "public", names(virginFunc), virginFunc)