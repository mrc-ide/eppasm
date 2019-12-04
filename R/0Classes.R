# Zero in file's name to make R loads this file first and make all classes
# available to use later on. Similar to `.hpp`. See the respective xClass for
# its methods

# EPP parameter super class: TODO: fp as an R6?
#' @importFrom methods new
eppFP <- R6::R6Class("eppfp", class=F, cloneable=F, portable=F, lock_objects=F,
    public = list(
        p = NULL,
        initialize = function(fp) {
            p <<- fp[-1]
            list2env(fp$ss, self)
        }
    )
)

# POP CLASS
# -----------------------------------------------------------------------------
# preferablly all fields should have a leading character
popEPP <- R6::R6Class("popepp", class=F, cloneable=F, portable=F, inherit=eppFP,
    lock_objects=F,
    public = list(
        data              = "array",
        VERSION           = "character",
        MODEL             = "numeric",
        MIX               = "logical",
        year              = "numeric",
        birth_age         = "numeric", birth_agrp        = "numeric",
        prev15to49        = "vector",  incid15to49       = "vector",
        prev              = 0,         incid             = 0,
        entrantprev       = "vector",  pregprevlag       = "vector",
        pregprev          = "vector",  # mistook for pregprevlag
        birthslag         = "array",   infections        = "array",
        hivdeaths         = "array",   natdeaths         = "array",
        adj_prob_age      = "array",   hivp_entrants_out = "array",
        incrate15to49_ts  = "vector",  prev15to49_ts     = "vector", # can be array
        data_db           = "array",    # virgin population
        data_active       = "array",
        artcov            = numeric(2), # initially no one on treatment
        prev_last         = 0, # last time step prevalence
        hiv_sx_prob       = "array",
        hiv_mr_prob       = "array",
        adj_prob_agr      = "array",
        rvec              = "vector",
        mr_prob_          = "array",
        p_hiv_death_age_  = "array",
        VIRGIN = "virEPP class",
        MSM    = "msmEPP class",
        FSW    = "fswEPP class",
        CLIENT = "cliEPP class"
        )
)

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