fit.mod <- function(obj, fplist, epp=FALSE, B0 = 1e5, B = 1e4, B.re = 3000, number_k = 500){
  ## fplist : updates to fixed parameters (fp) object to specify fitting options

  likdat <<- attr(obj, 'likdat')  # put in global environment for IMIS functions.
  if(epp)
    fp <<- update(attr(obj, 'eppfp'), list=fplist)
  else
    fp <<- update(attr(obj, 'specfp'), list=fplist)

  ## If IMIS fails, start again [warning: will produce infinite loop if there's a bug...]
  fit <- try(stop(""), TRUE)
  while(inherits(fit, "try-error")){
    start.time <- proc.time()
    fit <- try(IMIS(B0, B, B.re, number_k))
    fit.time <- proc.time() - start.time
  }

  fit$time <- fit.time
  fit$fp <- fp
  fit$likdat <- likdat

  rm(fp, likdat, pos=.GlobalEnv)

  return(fit)
}
