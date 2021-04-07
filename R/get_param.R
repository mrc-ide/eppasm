#' Get named model posterior model parameters as a data frame
#' 
get_param <- function(fit){

  val <- list()

  param <- lapply(seq_len(nrow(fit$resample)), function(i) fnCreateParam(fit$resample[i,], fit$fp))

  if(fit$fp$eppmod == "rlogistic"){
    val$rlog <- fit$resample[ , 1:4]
    colnames(val$rlog) <- c("log r0", "log r1", "log alpha", "tmid")
    val$iota_par <- fit$resample[ , 4 + 1]
  } else if(fit$fp$eppmod == "rhybrid"){
    val$rlog <- fit$resample[ , 1:4]
    colnames(val$rlog) <- c("log r0", "log r1", "log alpha", "tmid")
    val$rw <- fit$resample[ , 4 + 1:fit$fp$rt$n_rw]
    val$iota_par <- fit$resample[ , 4 + fit$fp$rt$n_rw + 1]
  } else if(fit$fp$eppmod == "rspline") {
    val$u <- fit$resample[ , 1:fit$fp$numKnots]
    val$beta <- t(sapply(param, "[[", "beta"))
  } else if(fit$fp$eppmod == "rtrend") {
    val$t0 <- vapply(param, "[[", numeric(1), "tsEpidemicStart")
    val$t1 <- vapply(lapply(param, "[[", "rtrend"), "[[", numeric(1), "tStabilize")
    val$r0 <- vapply(lapply(param, "[[", "rtrend"), "[[", numeric(1), "r0")
    val$beta <- t(vapply(lapply(param, "[[", "rtrend"), "[[", numeric(4), "beta"))
  }

  if(fit$fp$eppmod == "rtrend")
     val$iota <- rep(fit$fp$iota, length(param))
  else
    val$iota <- sapply(param, "[[", "iota")

  if("log_frr_adjust" %in% names(param[[1]]))
    val$log_frr_adjust <- sapply(param, "[[", "log_frr_adjust")
  
  if("v.infl" %in% names(param[[1]]))
    val$log_vinfl <- log(sapply(param, "[[", "v.infl"))

  if("ancbias" %in% names(param[[1]]))
    val$ancbias <- sapply(param, "[[", "ancbias")

  if("ancrtcens.vinfl" %in% names(param[[1]]))
    val$ancrtcens.vinfl <- sapply(param, "[[", "ancrtcens.vinfl")

  if("ancrtsite.beta" %in% names(param[[1]]))
    val$ancrtsite.beta <- sapply(param, "[[", "ancrtsite.beta")

  do.call(data.frame, val)
}
