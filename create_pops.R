
library(data.table)
library(readxl)
library(dplyr)
devtools::load_all("~/eppasm")

in.path = paste0("C:/Users/DeepaJahagirdar/Avenir Health Dropbox/Avenir Shared Drive/DataSets/UN POP/UN 2022 Revision/")
  
#All pops downloaded from: XX

#Read in pop
pop = read_xlsx(paste0(in.path,"ALL_MALE.xlsx")) %>% data.table()
pop = pop[12:nrow(pop)]
names = as.character(pop[1,])
colnames(pop) <- names
pop = pop[2:nrow(pop)]
setnames(pop,c("Region, subregion, country or area *","ISO3 Alpha-code","Type","Year"),
         c("country","iso3","type","year"))
ages = paste0(1:99)
pop = pop[,c("country","iso3","type","year",ages),with=FALSE]
pop = pop[type == "Country/Area"]
pop = pop[year %in% 1970:2019]
pop = melt(pop, id.vars = c("country","iso3","type","year"))
x1 = split(pop,c(pop$country))
x2 = lapply(x1,function(x) split(x, x$year))

unique(pop$year)
dim(subp[[1]]$N)
subp[[1]]$N[1,,]
