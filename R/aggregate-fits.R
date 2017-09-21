create_aggr_input <- function(inputlist){

  val <- inputlist[[1]]

  eppdlist <- lapply(inputlist, attr, "eppd")

  ## Only keep intersecting years
  anc.prev.list <- lapply(eppdlist, "[[", "anc.prev")
  anc.prev.list <- lapply(anc.prev.list, function(x) x[,Reduce(intersect, lapply(anc.prev.list, colnames))])

  anc.n.list <- lapply(eppdlist, "[[", "anc.n")
  anc.n.list <- lapply(anc.n.list, function(x) x[,Reduce(intersect, lapply(anc.n.list, colnames))])

  ## !! NEEDS ANC-RT DATA ADDED

  attr(val, "eppd") <- list(anc.used = do.call(c, lapply(eppdlist, "[[", "anc.used")),
                            anc.prev = do.call(rbind, anc.prev.list),
                            anc.n = do.call(rbind, anc.n.list))
  ## attr(val, "likdat") <- list(anclik.dat = with(attr(val, "eppd"), anclik::fnPrepareANCLikelihoodData(anc.prev, anc.n, anc.used, attr(val, "specfp")$ss$proj_start)))
  ## attr(val, "likdat")$lastdata.idx <- max(unlist(attr(val, "likdat")$anclik.dat$anc.idx.lst))
  ## attr(val, "likdat")$firstdata.idx <- min(unlist(attr(val, "likdat")$anclik.dat$anc.idx.lst))

  artnumperc <- !attr(inputlist[[1]], "specfp")$art15plus_isperc
  artnumlist <- lapply(lapply(inputlist, attr, "specfp"), "[[", "art15plus_num")

  art15plus_num <- artnumlist[[1]]
  art15plus_num[artnumperc] <- Reduce("+", lapply(artnumlist, "[", artnumperc))

  attr(val, "specfp")$art15plus_num <- art15plus_num

  return(val)
}


fnAddHHSLikDat <- function(obj){
  objcountry <- attr(obj, "country")
  fp <- attr(obj, "specfp")
  anchor.year <- as.integer(floor(min(fp$proj.steps)))

  attr(obj, "eppd")$hhs <- subset(prev.15to49.nat, country==objcountry)
  attr(obj, "eppd")$hhsage <- subset(prev.agesex.nat, country==objcountry)
  ## attr(obj, "eppd")$sibmx <- subset(sib.mx.tips, country==objcountry)

  return(obj)
}
