### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/eppasm-1/")
## Packages
library(data.table); library(mvtnorm); library(survey)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
	run.name <- args[1]
	loc <- args[2]
	proj.end <- as.integer(args[3])
	i <- as.integer(Sys.getenv("SGE_TASK_ID"))
} else {
	run.name <- "181126_test"
	loc <- "MWI"
	proj.end <- 2019
	i <- 1
}

### Arguments
start.year <- 1970
stop.year <- 2019
trans.params <- TRUE
gbd.pop <- TRUE
art.sub <- FALSE
prev.sub <- FALSE
collapse <- TRUE
anc.prior <- TRUE
no.anc <- FALSE
eq.prior <- TRUE
anc.backcast <- FALSE
num.knots <- 7
popadjust <- NULL
popupdate <- TRUE
use_ep5 = FALSE


### Paths
input.dir <- paste0()
out.dir <- paste0('/ihme/hiv/epp_output/gbd19/', run.name, "/", loc)
out.path <- paste0(out.dir, "/results", i, ".RData")
pdf.path <- paste0(out.dir, "/test_results", i, ".pdf")

### Functions
library(mortdb, lib = "/home/j/WORK/02_mortality/shared/r")
setwd(code.dir)
devtools::load_all()

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
## Prep data and collapse location subpopulations
dt <- prepare_spec_fit_gbd(loc, collapse, i, proj.end, gbd.pop, popadjust, popupdate, use_ep5)

## Substitute IHME data
# Prevalence surveys
if(prev.sub) {
	if((collapse & length(dt) == 1) | grepl("IND_", loc)) {
		print("Substituting prevalence surveys")
		dt <- sub.prev(loc, dt)	
	}
}

# # ANC data
# high.risk.list <- loc.table[epp == 1 & collapse_subpop == 0 & !grepl("ZAF", ihme_loc_id) & !grepl("KEN", ihme_loc_id), ihme_loc_id]
# ken.anc.path <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/data/prepped/kenya_anc_map.csv")
# ken.anc <- fread(ken.anc.path)
# no.anc.ken <- setdiff(loc.table[epp == 1 & grepl("KEN", ihme_loc_id), ihme_loc_id], ken.anc$ihme_loc_id)
# if(loc %in% c(high.risk.list, "PNG", no.anc.ken)) {
# 	anc.backcast <- F
# }
# if(anc.backcast) {
# 	dt <- sub.anc(loc, dt)
# }

# # Transition parameters
# if(trans.params) {
# 	dt <- extend.trans.params(dt, start.year, stop.year)
# 	dt <- sub.off.art(dt, loc, i)
# 	dt <- sub.on.art(dt, loc, i)
# 	dt <- sub.cd4.prog(dt, loc, i)
# }

## Fit model
fit <- list() 
for(subpop in names(dt)) {
	print(subpop)
	# if(subpop == "HSH") {
	# 	attr(dt[[subpop]], "likdat")$hhslik.dat[1:2, "used"] <- FALSE
	# }
	if(all(!attr(dt[[subpop]], "eppd")$anc.used)) {
        no.anc <<- TRUE
	} else {
		no.anc <<- FALSE
	}
	attr(dt[[subpop]], "eppd")$anc.used[1] <- FALSE
	attr(dt[[subpop]], "eppfp")$artelig.idx.ts <- as.integer(attr(dt[[subpop]], "eppfp")$artelig.idx.ts)
	# if(anc.prior) {
	# 	set.anc.prior(loc, subpop)
	# }
	fit[[subpop]] <- fitmod(dt[[subpop]], eppmod = 'rhybrid', rw_start = 2010,  B0=1e3, B=1e2, opt_iter=1:2*5, number_k = 50)
}


## Prepare output
result <- lapply(fit, simfit.gbd)
# result <- prep_epp_output(fit)
dir.create(out.dir, showWarnings = F, recursive = T)
save(result, file = out.path)

## Aggregate subpopulations to national and write prevalence and incidence draws
years <- unique(floor(result[[1]]$fp$proj.steps))
nat.data <- nat.draws(result)
var_names <- sapply(1:ncol(nat.data$prev), function(a) {paste0('draw',a)})
out_data <- lapply(nat.data, data.frame)
for (n in c('prev', 'incid', 'art', 'art_num', 'pop')) {  
  names(out_data[[n]]) <- var_names
  out_data[[n]]$year <- years
  col_idx <- grep("year", names(out_data[[n]]))
  out_data[[n]] <- out_data[[n]][, c(col_idx, (1:ncol(out_data[[n]]))[-col_idx])]
  write.csv(out_data[[n]], paste0(out.dir, "/results_", n , i, ".csv"), row.names=F)
}

## Plot results
plot.fit(result, pdf.path, nat.data)

### END