### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
root <- ifelse(windows,"J:/","/home/j/")
user <- ifelse(windows, Sys.getenv("USERNAME"), Sys.getenv("USER"))
code.dir <- paste0(ifelse(windows, "H:", paste0("/homes/", user)), "/HIV/")

## Packages
library(data.table); library(survey)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
	run.name <- args[1]
} else {
	run.name <- paste0(substr(gsub("-","",Sys.Date()),3,8), "_backcast2")
}

### Paths
out.dir <- paste0("/ihme/hiv/epp_output/gbd17/", run.name, "/")
dir.create(out.dir, showWarnings = F)
out.path <- paste0(out.dir, "prev_surveys.csv")
supp.survey.path <- paste0(root, "WORK/04_epi/01_database/02_data/hiv/04_models/gbd2015/02_inputs/supplement_survey_data_2017.csv")

### Functions
source(paste0(root, "Project/Mortality/shared/functions/get_locations.r"))

### Tables
loc.table <- data.table(get_locations(hiv_metadata = T))

### Code
## GBD locs
gbd.locs <- loc.table$ihme_loc_id

#bring in geospatial microdata data
geos_dir <- paste0(root, "LIMITED_USE/LU_GEOSPATIAL/geo_matched/hiv_gbd/")
versions <- grep("[[:digit:]]*\\.[[:digit:]]*\\.[[:digit:]]*", list.dirs(geos_dir, full.names=F), value=T)
newest <- versions[which.max(as.Date(versions, format="%m.%d.%y"))]
load(dir(paste0(geos_dir, newest), pattern=".Rdata", full.names=T)[1])
setnames(gbd_all, "country", "iso3")
data3 <- gbd_all[between(age_year, 15, 49) & !is.na(iso3) & !is.na(hiv_test) & !is.na(hiv_weight),]
# Kenya
data3[grepl("KEN", iso3) & !(admin_2_id == "KEN" | admin_2_id == "" | is.na(admin_2_id)), iso3 := admin_2_id]
data3 <- data3[iso3 != "KEN"]

# South Africa
data3[grepl("ZAF", iso3) & !is.na(admin_1_id), iso3 := admin_1_id]

# Ethiopia
data3[grepl("ETH", iso3) & !is.na(admin_1_id), iso3 := admin_1_id]

# India
data3[grepl("IND", iso3) & !is.na(admin_1_id), iso3 := admin_1_id]

# Repeat India states U/R
ind.copy <- copy(data3[grepl("IND", iso3)])
ind.copy[, iso3 := admin_1_urban_id]
data3 <- rbind(data3, ind.copy)
data3[,loc_year := paste0(iso3,"_",year)] 

## Nigeria
# data3[grepl('NGA', iso3) & !is.na(admin_1_id), iso3 := admin_1_id]


#bring in report data (not extracted yet)
supp.survey <- fread(supp.survey.path)[iso3 %in% gbd.locs & outlier == 0]
supp.survey[, outlier := NULL]

## process microdata
#need this option unfortunately
options(survey.lonely.psu="adjust")
data4 <- rbindlist(
	lapply( unique(data3$loc_year),
		function(loc.year){
			print(loc.year)
			data5 <- data3[loc_year == loc.year]
			loc <- unique(data5$iso3)
			if(all(data5$hiv_test == 0)) {
				# impute a half positive observation to get uncertainty for 0 prevalence
				hold <- copy(data5[1,])
				hold$hiv_test <- 0.5
				data5 <- rbind(data5, hold)
			}
			if(nrow(data5[is.na(strata),]) > 0){
				print(paste0((nrow(data5[is.na(strata),])/nrow(data5))*100,"% missing strata for ",loc.year))
				if(nrow(data5[!is.na(strata),]) == 0) {
					data6 <- data5
					s <- svydesign(ids = ~psu, weights = ~hiv_weight,data = data6)
				} else {
					data6 <- data5[!is.na(strata),]
					s <- svydesign(ids = ~psu, strata = ~strata, weights = ~hiv_weight,data = data6,check.strata = TRUE)	                
				}
				t <- svymean(~hiv_test,s)
				data5 <- data5[,list(year = median(int_year), prev = weighted.mean(hiv_test, hiv_weight, na.rm = T), n = .N),
				by='iso3,nid,survey_name']
				d <- data.table(iso3 = loc, year = data5[,"year", with = F][[1]],prev =data5[,"prev", with = F][[1]], se = SE(t)[[1]], n = nrow(data6)  )
				return(d)
			} else { 
				if(length(unique(data5$psu)) == 1) {
					data6 <- data5
					s <- svydesign(id = ~hh_id, strata = ~strata, weights = ~hiv_weight,data = data6)
				} else {
					s <- svydesign(ids = ~psu, strata = ~strata, weights = ~hiv_weight, data = data5, check.strata = TRUE)
				}	          	
				t <- svymean(~hiv_test,s)
				data13 <- data5[,list(year = median(int_year), prev=weighted.mean(hiv_test, hiv_weight, na.rm=T), n=.N),
				by='iso3,nid,survey_name']
				d <- data.table(iso3 = loc, year = data13[,"year", with = F][[1]],prev =coef(t)[[1]], se = SE(t)[[1]], n = nrow(data5)  )
				return(d)
			}

		}
	)
)
nid.dt <- unique(data3[, .(int_year, iso3, nid)])
nid.dt[, year := as.numeric(int_year)]
temp <- merge(data4, nid.dt, by = c("year", "iso3"), all.x = T)
out.dt <- rbind(temp, supp.survey, fill = T)

# Check for duplicates
loc.years <- data.table(melt(table(out.dt[, .(iso3, year)])))
loc.years[value == 2]

# Outlier surveys
out.dt <- out.dt[!(grepl("ZAF", iso3) & year %in% c(2002, 2008))]
out.dt <- out.dt[!(iso3 %in% c("KEN_35623", "KEN_35662") & year == 2008)]

write.csv(out.dt, out.path, row.names = F)


### End