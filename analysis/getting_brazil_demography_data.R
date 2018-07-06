########################################################################################################################
## Sorting out population file for input into the model ################################################################
########################################################################################################################

total_pop <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/WPP2017_TotalPopulationBySex.csv")

brazil_pop <- total_pop[total_pop$Location == "Brazil",]

brazil_pop_period <- subset(brazil_pop, brazil_pop$Time > 1969 & brazil_pop$Time <= 2021)

write.csv(brazil_pop_period,
          file = "C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_WPP_brazil_data.csv",row.names = F)
## Cheeky pop plot 

pop_estimates <- subset(brazil_pop_period, brazil_pop_period$Variant == "Low" | brazil_pop_period$Variant == "Medium" |
                          brazil_pop_period$Variant == "High" )
a<-ggplot(data = pop_estimates,aes(x=Time,y=PopTotal,group=Variant)) + geom_line(aes(colour=Variant),size=1.05)

plot(pop_estimates[pop_estimates$Variant=="Low",])

#######################################################################################################################
## Now lets get the total pop in there ################################################################################
#######################################################################################################################

pop_by_age <- read.csv("c:/Users/josh/Dropbox/hiv_project/population_data/WPP2017_PopulationBySingleAgeSex.csv")

brazil_pop_by_age <- subset(pop_by_age,pop_by_age$Location == "Brazil")
brazil_pop_by_age

brazil_pop_age_period <- subset(brazil_pop_by_age, brazil_pop_by_age$Time >= 1970 & brazil_pop_by_age$Time <= 2021)

ages <- as.character(15:80)
ages[67] <- "80+"

brazil_adult_pop <- NULL
for(i in 1:nrow(brazil_pop_age_period)){
  row_current <- brazil_pop_age_period[i,]
  if(row_current$AgeGrp %in% ages){
    brazil_adult_pop <- rbind(brazil_adult_pop,row_current)
  }
}

write.csv(brazil_adult_pop,
          file = "C:/Users/josh/Dropbox/hiv_project/population_data/trimed_pop_by_age_ADULT_BRAZIL.csv",row.names = F)

#######################################################################################################################
## Lets work out the paed pop per year ################################################################################
#######################################################################################################################

paed_pop <- as.character(0:14)

brazil_pop_age_period
brazil_paed_df <- NULL

for(i in 1970:2021){
  df <- brazil_pop_age_period[brazil_pop_age_period$Time == i, ]
  ages_want_df <- subset(df,df$AgeGrp == paed_pop)
  
  paed_pop_row <- cbind.data.frame(ages_want_df[1,1:9],sum(ages_want_df[,10]),sum(ages_want_df[,11]),sum(ages_want_df[,12]))
  
  names(paed_pop_row) <- names(brazil_pop_age_period)
  
  brazil_paed_df <- rbind.data.frame(brazil_paed_df,paed_pop_row)
  
  }
brazil_paed_df$AgeGrp <- rep("0-14",nrow(brazil_paed_df))
brazil_paed_df$AgeGrpSpan <- rep("15"
                                 ,nrow(brazil_paed_df))


write.csv(brazil_paed_df,
          file = "C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_paed_data_grouped_BRAZIL.csv")

#####################################################################################################################
## Now we will get the survival probability data ####################################################################
#####################################################################################################################

rm(pop_by_age)
rm(total_pop)

sx_data <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/WPP2017_LifeTable.csv")

sx_brazil_data <- subset(sx_data, Location =="Brazil")

sx_brazil_data


#####################################################################################################################
## Deaths data ######################################################################################################
#####################################################################################################################

AIDS_deaths <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/Total_deaths_2017_2.csv")
lines(cl_csavrd$aids_deaths)
plot(cl_csavrd$aids_deaths,type="l",ylim=c(0,550))
lines(nl_csavrd$aids_deaths,type="l",col="blue")
abline(v=25,col="red")
#####################################################################################################################
## NOw lets get the births merely from the numbers aged 0 in a year .... ############################################
#####################################################################################################################

tot_pop <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/WPP2017_PopulationBySingleAgeSex.csv")

brazil_pop <-subset(tot_pop, Location == "Brazil")

brazil_period <- subset(brazil_pop, Time > 1969 & Time <2022)

births_brazil <- subset(brazil_period, AgeGrp == "0")

write.csv(births_brazil,
          file="C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_birth_data_from_age_0_total_pop_BRAZIL.csv")

#####################################################################################################################
## Now lets get the entrant pop, i.e those aged 14 in the previous year #############################################
#####################################################################################################################

brazil_14_period <- subset(brazil_pop, Time > 1968 & Time < 2021)

entrant_pop_brazil <- subset(brazil_14_period, AgeGrp == "14")

entrant_pop_brazil_absolute <- entrant_pop_brazil$PopTotal

write.csv(entrant_pop_brazil,
          file = "C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_entrant_pop_from_14_years_BRAZIL.csv")

#####################################################################################################################
## Now lets try and get the fp object with the brazillian data ######################################################
#####################################################################################################################

brazil_deaths <- read.csv("C:/Users/josh/Dropbox/hiv_project/brazil_mortality_data/Total_deaths_2017_2.csv")

brazil_births <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_birth_data_from_age_0_total_pop_BRAZIL.csv")

brazil_entrant <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/trimmed_entrant_pop_from_14_years_BRAZIL.csv")

brazil_pop <- read.csv("C:/Users/josh/Dropbox/hiv_project/population_data/trimed_pop_by_age_ADULT_BRAZIL.csv")

target_pop <- array(0,c(66,2,52))

for(i in 1970:2021){
  year_pop <- brazil_pop[brazil_pop$Time==i,]
  
  slice <- i - 1969
  
  target_pop[,1,slice] <- year_pop$PopMale * 1000
  target_pop[,2,slice] <- year_pop$PopFemale * 1000
  
  
}

births <- brazil_births$PopTotal

names(births) <- c(1970:2021)

births <- births * 1000

basepop <- array(0,c(66,2))

basepop[,1] <- brazil_pop[brazil_pop$Time == 1970, 10]

basepop[,2] <- brazil_pop[brazil_pop$Time == 1970, 11]



cl_pjnz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Chile_2017_final.pjnz"

cl_fp <- prepare_directincid(cl_pjnz)

brazil_fp <- cl_fp

brazil_fp$basepop <- basepop

brazil_fp$births <- births

brazil_fp$targetpop <- target_pop

brazil_pnjz <- "C:/Users/josh/Dropbox/hiv_project/jeff_eppasm_data/Brazil_2017_final.PJNZ"

brazil_fp <- prepare_directincid(brazil_pnjz)
brazil_fp$t_diagn_start <- 16L


brazil_simmod <- simmod.specfp(brazil_fp)
