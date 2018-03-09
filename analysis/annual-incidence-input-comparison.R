## devtools::install_github("mrc-ide/eppasm@directincid", auth_token="<YOUR TOKEN>")

devtools::load_all("~/Documents/Code/R/eppasm-csavr/") # @directincid

pjnz <- c(list.files("~/Documents/Data/Spectrum files/2017 final/WCENA", "PJNZ$", full.names=TRUE),
          list.files("~/Documents/Data/Spectrum files/2017 final/LA", "PJNZ$", full.names=TRUE),
          list.files("~/Documents/Data/Spectrum files/2017 final/CAR", "PJNZ$", full.names=TRUE))

## Prepare model inputs with annual incidence input
inputs <- lapply(pjnz, prepare_directincid)
names(inputs) <- sapply(inputs, attr, "country")

## Simulate the model
mod <- lapply(inputs, simmod)

## Read Spectrum file ouputs for comparison
specres <- lapply(pjnz, read_hivproj_output)
names(specres) <- sapply(specres, attr, "country")

## Plot comparison

pdf("~/Documents/Code/R/eppasm-csavr/dev/wcena-la-car-comparison_2017-10-30.pdf", h=11, w=8.5, pointsize=9)
par(mfrow=c(4,3), mgp=c(2, 0.5, 0), tcl=-0.25, mar=c(2, 3, 2.5, 1), cex=1.0)
####n
for(ii in seq_along(mod)){
  country <- names(mod)[ii]
  modi <- mod[[ii]]
  spri <- specres[[ii]]
####
  matplot(as.integer(dimnames(spri$totpop)[[3]]),
          cbind(colSums(attr(modi, "infections"),,2),
                colSums((spri$newinf.m + spri$newinf.f)[4:17,])),
          type="l", col=c("red3", "grey30"), lwd=2, lty=1:2,
          xlim=c(1980, 2021), xlab="", ylab="new infections",
          main=paste0(country, "; new infections, 15+"))
  legend("bottomright", c("EPPASM", "UNAIDS17"), col=c("red3", "grey30"),
         lwd=2, lty=1:2, cex=0.8)
  ##
  matplot(as.integer(dimnames(spri$totpop)[[3]]),
          cbind(colSums(modi[,,2,],,2), colSums(spri$hivpop[16:81,,],,2)),
          type="l", col=c("red3", "grey30"), lwd=2, lty=1:2,
          xlim=c(1980, 2021), xlab="", ylab="HIV population",
          main=paste0(country, "; PLHIV, 15+"))
  legend("bottomright", c("EPPASM", "UNAIDS17"), col=c("red3", "grey30"),
         lwd=2, lty=1:2, cex=0.8)
  ##
  matplot(as.integer(dimnames(spri$totpop)[[3]]),
          cbind(colSums(attr(modi, "hivdeaths"),,2), colSums(spri$hivdeaths[16:81,,],,2)),
          type="l", col=c("red3", "grey30"), lwd=2, lty=1:2,
          xlim=c(1980, 2021), xlab="", ylab="AIDS deaths",
          main=paste0(country, "; AIDS deaths, 15+"))
  legend("bottomright", c("EPPASM", "UNAIDS17"), col=c("red3", "grey30"),
         lwd=2, lty=1:2, cex=0.8)
}
####
dev.off()


## Time to simulate the model for each Spectrum file (64 files)
library(microbenchmark)
microbenchmark(lapply(inputs, simmod))
