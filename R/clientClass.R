# client CLASS methods, see def in 0Classes.R
# -----------------------------------------------------------------------------
cliFunc <- c(
aging = function(ag_prob) {
    data[-c(1,pAG),,year] <<- data[-(pAG-1:0),,year-1]
    data[pAG,,year]       <<- data[pAG,,year-1] + data[pAG-1,,year-1] # open 
},

initiation = function(n_active) {
    data[,,year] <<- n_active[,m.idx,] * p$cli_init_rate[,,year] # by age/H+/-
},

turn_over = function() {
    data[,,year] <<- data[,,year] * (1 - p$cli_turn_over[,,year]) # by age/H+/-
},

deaths = function() {
    data[,,year] <<- sweep(data[,,year], 1, p$Sx[,m.idx,year], '*')
},

migration = function(migration_pr) {
    data[,,year] <<- sweep(data[,,year], 1, migration_pr[,m.idx], '*')
},

adjust_pop = function(adj_prob) {
    data[,,year] <<- sweep(data[,,year], 1, adj_prob[,m.idx,year], '*')
}, 

n_all = function(dt=0, only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    sum(data[age_range,,year]) -
        sum(data[age_range[1],,year]) * dt +
            sum(data[tail(age_range)+1,,year]) * dt
},

n_HIV = function(dt=0, only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    sum(data[age_range, hivp.idx, year]) -
        sum(data[age_range[1], hivp.idx, year]) * dt +
            sum(data[tail(age_range)+1, hivp.idx, year]) * dt
},

n_NEG = function(only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    data[age_range, hivn.idx, year]
},

n_IDU = function(only_1549=TRUE) {
    age_range = if (only_1549) p.age15to49.idx else 1:pAG
    data[age_range,,year] * p$idu_cli_pro[age_range,]
}
)

setMembers(cliEPP, "public", names(cliFunc), cliFunc)