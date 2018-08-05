library(devtools)
system.file("extdata")

#' This script adds some PNJZ files to the package for software testing and illustration purposes.
#' Files are taken from UNAIDS estimates files produced by national HIV estimates teams for UNAIDS.
#' However, files may not include the most current data and have been modified to test or demonstrate
#' specific software features and thus **should NOT be taken as estimates for given countries.**
#'
#' To access data and estimates for these and other countries, please visit http://aidsinfo.unaids.org/
#' or http://www.unaids.org/en/dataanalysis/datatools/spectrum-epp to request the most recent Spectrum
#' estimates files.


#' ## Botswana 2017 file
file.copy("~/Documents/Data/Spectrum files/2017 final/2017 Final Spectrum files/SSA/Botswana_29_05_2017 updated BF.PJNZ",
          "~/Documents/Code/R/eppasm/inst/extdata/testpjnz/Botswana2017.PJNZ")

fp <- "~/Documents/Code/R/eppasm/inst/extdata/testpjnz/Botswana2017.PJNZ"
unzip(fp, list=TRUE)

## Remove files to reduce size
zip(fp, grep("bm2$", unzip(fp, list=TRUE)$Name, value=TRUE), flags="-d")
zip(fp, grep("SPU$", unzip(fp, list=TRUE)$Name, value=TRUE), flags="-d")
zip(fp, grep("UAC$", unzip(fp, list=TRUE)$Name, value=TRUE), flags="-d")


#' ## Dominican Republic 2017 file

file <- "~/Documents/Code/R/eppasm/inst/extdata/testpjnz/DominicanRepublic2017.PJNZ"

file.copy("~/Documents/Data/Spectrum files/2017 final/CAR/Dominican Republic_2017_final.PJNZ", file)

unzip(fp, list=TRUE)

## Remove files to reduce size
zip(file, grep("bm2$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")
zip(file, grep("SPU$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")
zip(file, grep("UAC$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")


#' ## Dominican Republic 2017 file

file <- "~/Documents/Code/R/eppasm/inst/extdata/testpjnz/Netherlands2017.PJNZ"

file.copy("~/Documents/Data/Spectrum files/2017 final/WCENA/Netherlands_2017_final.PJNZ", file)

unzip(fp, list=TRUE)

## Remove files to reduce size
zip(file, grep("bm2$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")
zip(file, grep("SPU$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")
zip(file, grep("UAC$", unzip(file, list=TRUE)$Name, value=TRUE), flags="-d")
