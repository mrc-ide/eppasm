#' Read Spectrum _pop1.xlsx export file
#'
#' Spectrum has a debug mode option to output a file recording
#' the entire state space when a model simulation is computed.
#' This function parse the *_pop1.xlsx file into long-format
#' data.frame.
#'
#' @param pop1file path to *_pop1.xlsx export file
#' @param country string giving the identifying country / region
#'    name to be added as a column in the output file
#' @param years integer vector of years to extract
#'
#' @return A data.frame consisting of pop1 file in long format
#'   
#' @details
#'
#' TODO: tidy up the dependencies in this function: reshape, pryr, readxl.
#' 
#' TODO: Insert instructions for how to enter Spectrum debug mode.
#' 
#' `years` must be in the file as sheets, but this is not currently checked.
#' It will fail hard. Inserting a check for this would make it fail softer.
#'
#' @export
#' 
read_pop1 <- function(pop1file, country, years = 2000:2021){

  print(country)

  ## Read population out of Spectrum file
  pop_list <- list()
  for(yr in years){
    pop1 <- as.data.frame(readxl::read_excel(pop1file, sheet=as.character(yr)))
    pop1 <- pop1[-1,]
    names(pop1)[1:3] <- c("sex", "cd4", "artdur")
    pop1 <- reshape(pop1, idvar=c("sex", "cd4", "artdur"), varying=c(0:79, "80+"), times=c(0:79, "80+"), timevar="age", v.names="pop", direction="long")
    pop1$year <- as.integer(yr)
    pop1$age <- as.integer(gsub("\\+", "", pop1$age))
    pop1$cd4 <- as.integer(gsub("CD4 Count: ", "", pop1$cd4))
    pop1$artdur <- as.integer(gsub("Duration: ", "", pop1$artdur))
    pop_list[[as.character(yr)]] <- pop1
  }
  pop <- plyr::rbind.fill(pop_list)
  pop$sex <- gsub("ales", "ale", pop$sex)
  pop <- subset(pop, cd4 != 0 & artdur != 0 & pop > 0)

  pop$country <- country

  return(pop)
}
