obj = dt
epp=FALSE
B0 = 1e5
B = 1e4
number_k = 500
opt_iter=0, 
sample_prior=eppasm:::sample.prior
prior=eppasm:::prior
likelihood=eppasm:::likelihood
fp <- update(attr(obj, 'specfp'))
fp$eppmod = epp.mod  
eppd <- attr(obj, "eppd")

has_ancrtsite <- exists("ancsitedat", eppd) && any(eppd$ancsitedat$type == "ancrt")
has_ancrtcens <- !is.null(eppd$ancrtcens) && nrow(eppd$ancrtcens)
fp$ancrt <- "site"
likdat <- prepare_likdat(eppd, fp)
fp$ancsitedata <- as.logical(nrow(likdat$ancsite.dat$df))
fp$SIM_YEARS <- as.integer(max(likdat$ancsite.dat$df$yidx,
                               likdat$hhs.dat$yidx,
                               likdat$ancrtcens.dat$yidx,
                               likdat$hhsincid.dat$idx,
                               likdat$sibmx.dat$idx,
                               max(as.integer(colnames(likdat$vr))) - min(as.integer(colnames(likdat$vr))) + 1))

fp$proj.steps <- seq(fp$ss$proj_start+0.5, fp$ss$proj_start-1+fp$SIM_YEARS+0.5, by=1/fp$ss$hiv_steps_per_year)
tsEpidemicStart <- if(epp) fp$tsEpidemicStart else fp$ss$time_epi_start+0.5
fp <- prepare_rhybrid(fp)
fp$logitiota <- TRUE

fp$incidmod <- "eppspectrum"

fp <- prepare_irr_model(fp)
X_k <- sample_prior(B0, fp)
cov_prior = cov(X_k)        # estimate of the prior covariance

## Locations and covariance of mixture components
center_all <- list()
sigma_all <- list()

ll_all <- numeric()
lprior_all <- numeric()
mix_all <- numeric()
n_k <- integer()

X_all <- matrix(0, B0+B*number_k, ncol(X_k))
n_all <- 0
stat <- matrix(NA, number_k, 7)

idx_exclude <- integer()  ## inputs to exclude from initial points for optimization step
k = 1
x = fread('/share/hiv/epp_output/gbd19/190205_nobackcast_agesexdat/MWI/theta_99.csv')
theta = x$theta
theta.last <<- theta
fp <- update(fp, list=fnCreateParam(theta, fp))
