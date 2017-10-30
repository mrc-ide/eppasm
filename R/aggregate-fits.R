create_aggr_input <- function(inputlist){

  val <- inputlist[[1]]

  eppdlist <- lapply(inputlist, attr, "eppd")

  ## Only keep intersecting years
  anc.prev.list <- lapply(eppdlist, "[[", "anc.prev")
  anc.prev.list <- lapply(anc.prev.list, function(x) x[,Reduce(intersect, lapply(anc.prev.list, colnames))])

  anc.n.list <- lapply(eppdlist, "[[", "anc.n")
  anc.n.list <- lapply(anc.n.list, function(x) x[,Reduce(intersect, lapply(anc.n.list, colnames))])

  ancrtsite.prev.list <- lapply(eppdlist, "[[", "ancrtsite.prev")
  ancrtsite.prev.list <- lapply(ancrtsite.prev.list, function(x) x[,Reduce(intersect, lapply(ancrtsite.prev.list, colnames))])

  ancrtsite.n.list <- lapply(eppdlist, "[[", "ancrtsite.n")
  ancrtsite.n.list <- lapply(ancrtsite.n.list, function(x) x[,Reduce(intersect, lapply(ancrtsite.n.list, colnames))])


  ## aggregate census data across regions
  ancrtcens <- do.call(rbind, lapply(eppdlist, "[[", "ancrtcens"))
  if(!is.null(ancrtcens) && nrow(ancrtcens)){
    ancrtcens$x <- ancrtcens$prev * ancrtcens$n
    ancrtcens <- aggregate(cbind(x,n) ~ year, ancrtcens, sum)
    ancrtcens$prev <- ancrtcens$x / ancrtcens$n
    ancrtcens <- ancrtcens[c("year", "prev", "n")]
  }

  attr(val, "eppd") <- list(anc.used = do.call(c, lapply(eppdlist, "[[", "anc.used")),
                            anc.prev = do.call(rbind, anc.prev.list),
                            anc.n = do.call(rbind, anc.n.list),
                            ancrtsite.prev = do.call(rbind, ancrtsite.prev.list),
                            ancrtsite.n = do.call(rbind, ancrtsite.n.list),
                            ancrtcens = ancrtcens)

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
