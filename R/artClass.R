#' Set a list of methods to a R6 object
#'
#' @param a the object
#' @param where public, private or active?
#' @param x a named list,
#' @param y the corresponding functions to x
setMembers <- function(a, where="public", x, y) {
    if (length(x) > 1)
        mapply(function(x, y) a$set(where, x, y, overwrite=T), x, y)
    else 
        a$set(where, x, y, overwrite=T)
    invisible()
}

# ART GRAD
# -----------------------------------------------------------------------------
artEPP <- R6::R6Class("artepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F,
    public = list(
        year       = 1,
        MODEL      = "integer",
        data       = "array",
        data_db    = "array",
        gradART    = "array",
        gradART_db = "array",
        f_death     = "array",
        f_death_db  = "array",
        initialize = function(fp, MODEL) {
            super$initialize(fp)
            MODEL <<- MODEL
            data <<- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))            
            gradART <<- array(0, c(hTS, hDS, hAG, NG))
            f_death  <<- array(0, c(hTS, hDS, hAG, NG))
            if (MODEL==2) {
                data_db <<- data
                gradART_db <<- gradART
                f_death_db  <<- array(0, c(hTS, hDS, hAG, NG))
            }
        }
    )
)

artFns <- c(
aging = function(ag_prob) {
    data[,,,,year] <<- data[,,,,year-1]
    nARTup <- sweep(data[,,-hAG,,year], 3:4, ag_prob[-hAG,], "*")
    data[,,-hAG,,year] <<- data[,,-hAG,,year] - nARTup
    data[,,  -1,,year] <<- data[,,  -1,,year] + nARTup
    if (MODEL==2) {
        data_db[,,,,year] <<- data_db  [,,,,year-1]
        nARTup <- sweep(data_db[,,-hAG,,year], 3:4, ag_prob[-hAG,], "*")
        data_db[,,-hAG,,year] <<- data_db  [,,-hAG,,year] - nARTup
        data_db[,,  -1,,year] <<- data_db  [,,  -1,,year] + nARTup        
    }
},

add_entrants = function(artYesNo) {
    n_in <- sweep(p$paedsurv_artcd4dist[,,,year], 3, head(artYesNo, 2), "*")
    if (MODEL==1)
        data[,,1,,year] <<- data[,,1,,year] + n_in
    if (MODEL==2) # add to virgin then debut
        data_db[,,1,,year] <<- data_db[,,1,,year] + n_in
},

sexual_debut = function() {
    debut_art <- sweep(data_db[,,db_agr,,year], 3:4, p$db_pr, '*')
       data[,,db_agr,,year] <<-    data[,,db_agr,,year] + debut_art
    data_db[,,db_agr,,year] <<- data_db[,,db_agr,,year] - debut_art
},

deaths = function(survival_pr) {
    data[,,,,year] <<- sweep(data[,,,,year], 3:4, survival_pr, "*")
    if (MODEL==2)
        data_db[,,,,year] <<- sweep(data_db[,,,,year], 3:4, survival_pr, "*")
},

migration = function(migration_pr) {
    data[,,,,year] <<- sweep(data[,,,,year], 3:4, migration_pr, "*")
    if (MODEL==2)
        data_db[,,,,year] <<- sweep(data_db[,,,,year], 3:4, migration_pr, "*")
},

count_death = function() {
    f_death <<- p$art_mort * p$artmx_timerr[,year] * data[,,,,year]
    if (MODEL==2)
        f_death_db <<- p$art_mort * p$artmx_timerr[,year] * data_db[,,,,year]
},

grad_progress = function() {
    gradART[,,,] <<- 0 # reset every time step the gradient
    ## progression and mortality (HARD CODED 6 months duration)
    art_up <- 2.0 * data[1:2,,,,year]
    gradART[1:2,,,] <<- gradART[1:2,,,] - art_up
    gradART[2:3,,,] <<- gradART[2:3,,,] + art_up  
    # ART mortality
    gradART <<- gradART - f_death
    if (MODEL==2) {
        gradART_db[,,,] <<- 0 # reset every time step the gradient
        art_up <- 2.0 * data_db[1:2,,,,year]
        gradART_db[1:2,,,] <<- gradART_db[1:2,,,] - art_up
        gradART_db[2:3,,,] <<- gradART_db[2:3,,,] + art_up
        # ART mortality
        gradART_db <<- gradART_db - f_death_db
    }
},

art_dropout = function(hivpop) {
    now          <- data[,,,,year]
    hivpop$grad  <- hivpop$grad + colSums(now) * p$art_dropout[year]
    gradART     <<- gradART             - now  * p$art_dropout[year]
    if (MODEL==2) {
        now            <- data_db[,,,,year]
        hivpop$grad_db <- hivpop$grad_db + colSums(now) * p$art_dropout[year]
        gradART_db     <<- gradART_db            - now  * p$art_dropout[year]
    }
},

current_on_art = function() {
    art_curr <- colSums(data[,,,,year],,3) + colSums(gradART,,3) * DT
    if (MODEL==2) 
        art_curr <- art_curr +
                    colSums(data_db[,,,,year],,3) + colSums(gradART_db,,3) * DT
    art_curr
},

grad_init = function(artinit) {
    gradART[1,,,] <<- gradART[1,,,] + artinit / DT
    data[,,,,year] <<- data[,,,,year] +  DT * gradART
},

grad_db_init = function(artinit_db) {
    gradART_db[1,,,] <<- gradART_db[1,,,] + artinit_db / DT
    data_db[,,,,year] <<- data_db[,,,,year] +  DT * gradART_db
},

adjust_pop = function(adj_prob) {
    data[,,,,year] <<- sweep(data[,,,,year], 3:4, adj_prob, "*")
    if (MODEL==2)
        data_db[,,,,year] <<- sweep(data_db[,,,,year], 3:4, adj_prob, "*")
},

set_data = function(FUN="+", x, TS=T, DS=T, AG=T, NG=T, YEAR=NULL) {
    FUN <- match.fun(FUN)
    if (is.null(YEAR))
        data[TS, DS, AG, NG] <<- FUN(data[TS, DS, AG, NG], x)
    else
        data[TS, DS, AG, NG, YEAR] <<- FUN(data[TS, DS, AG, NG, YEAR], x)
}, 

get = function(YEAR=NULL, TS=T, DS=T, AG=T, NG=T) {
    'get does not change anything, just return the data, use set to change'
    if (is.null(YEAR))
        data[TS, DS, AG, NG, TRUE]
    else
        data[TS, DS, AG, NG, YEAR]
},

sweep_sex = function(FUN="*", x, year) {
    data[,,,,year] <<- sweep(data[,,,,year], 3:4, x, FUN)
})

setMembers(artEPP, "public", names(artFns), artFns)