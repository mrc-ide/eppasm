## eppasm 0.5.10

* Patch `read_epp_t0()` to parse XML file projection sets using same code as `epp:read_epp_data()`. This resolves parsing issue for Niger, Senegal, and Dominican Republic. Unsure why the XML is formed different for these cases.


## eppasm 0.5.9

* Update `read_hivproj_param()` to read incidence sex ratio from Spectrum v5.88, v5.89, and v5.90.

## eppasm 0.5.8

- Update read_incid_input() to use tag "<IncidenceByFit MV4>" in Spectrum DP file.
- Update read_hivproj_param() to parse incidence rate ratio fom new "<SexRatioByEpidPatt MV>" tag.
- Add model output `artinit` for number of ART initiations by CD4 category.
- Add model outputs `aidsdeaths_noart` and `aidsdeaths_noart` for number of AIDS deaths by 
  coarse age groups, stage of infection, and ART duration.
- Make ART stage durations (0-5 months, 6-11 months, 1+ years) input parameters `fp$ss$h_art_stage_dur` instead of hard coded values.
- Add function `spec_add_dimnames()` to add dimnames to simmod output arrays.

## eppasm 0.5.6

- Re-implement r-hybrid model and incorporate smooth transition from logistic function to random walk.
- Adjustable knot spacing for random walk component of r-hybrid model.
- Implement Spectrum options for new ART patient allocation.
- Implement option to scale mortality among untreated population by ART coverage (`fp$scale_cd4_mort = 1`)
- Separate time time varying ART mortality trends for ART duration <12 months and > 12 months.


## eppasm 0.5.5

- Allow specification equilibrium prior standard deviation as a fixed prior parameter.
- Revise logistic RW model to use fixed standard deviation.
- Rename rlogistic_rw -> rhybrid
- Add ANC site-level random effect estimates to tidy_output().
- Add get_pointwise_ll() to calculate pointwise log-likelihood for data or new data inputs.
- Add tidy_aggr() to generate summary outputs for aggregation of subnational fits.


## eppasm 0.5.4

- Updates for Spectrum v5.72:
- Update fertility rate ratio input parameters to allow ART FRR by age.
- Add time-varying mortality on ART


## eppasm 0.5.3

- More constrained  prior for log(iota) ~ Unif(log(1e-13), log(0.0025)), matching EPP software.
- Add posterior predictive outputs for site-level ANC data
- Incorporate an effect of percentage circumcised on male HIV incidence rate.


## eppasm 0.5.2

- Use Beer's coefficients to graduate 5-year age IRRs to single year.
- Add vignette to develop EPP-ASM to Spectrum workflow.


## eppasm 0.5.1

- Add tidy_output() function to generate long format of key output indicators.


## eppasm 0.5.0
- Implement age-specific prevalence among pregnant women.
- Refactor age-specific prevalence function.
- Add age-specific ART coverage function.

# eppasm 0.4.0
- Add annual direct incidence input option.

# eppasm 0.3.4
- Merge updates from v0.1.2 into v0.3.3

# eppasm 0.3.3

- Modularize code for calculating new infections.

- Implement age 15 entrant prevalence and ART coverage inputs.
- Implement logit transformation for sampling from uniform logiota prior distribution.
- Set age IRR penalty variance at a fixed constant value (9003)

- Implement piecewise-linear time-varying M:F incidence/transmission ratio
  and 15-24 incidence rate ratio. Knots fixed at 2002, 2007, and 2012.
- Modularize code for incidence rate ratio parameterization, now in infections.R.

- Updated default IRR priors to informative priors based on fits countries with several surveys.
- Add function to generate summary model outputs


# eppasm 0.3.2

- Integrate random effects variance parameter out of likelihood.
- Implement hybrid r-spline and RW model for r(t).
- Implement logistic model for log r(t) (rlogistic), and version with random walk (rlogistic_rw).
- Add optimization option to fitmod().

- Add option to fit without site-level ANC data.

# eppasm 0.3.0

- Rename package `eppasm`.
- Add IMIS function to eppasm package. In this version of IMIS, at each iteration a new mixture component is constructed either centered on input with greatest weight or based on optimizer. Allows optimizers to be run at arbitrary iterations, but doesn't implment multiple optimizers which might still be usful. Covariance for optimizer mixture component when Hessian is degenerate still needs further work.

- Change from dependency to imports epp package.
- Add ANC-RT data to aggregate national fits (prepare_national_fit() and create_aggr_input()).

- Implement dsamp() to calculate density for initial IMIS samples. This allows different initial sampling density from prior distribution.

- Specify knot locations for o-spline model
- Estimate HIV+:HIV- FRR adjustment.
- Specify knot locations for o-spline model

- Add functions to read country and subnational region from Spectrum PJN file.


# eppspectrum 0.1.6

- Bug fix: add check in ll_hhage_binom() that ldbinom() does not return NaN.
- Bug fix: in create_spectrum_fixpar(), handle case where ART initiation starts in different years for men and women.


# eppspectrum 0.1.5

- Implemented O-spline model (penalized B-spline) using mgcv to construct
- Added `r0logiotaratio` option to estimate ratio of r(t0) and log of initial seed, reduces correation. Not set as default.


# eppspectrum 0.1.3

- EPP model setup occurs in fitmod()
- Allow order of spline model penalty to vary (fp$rtpenord)
- Allow number of knots to vary (fp$numKnots)
- Return spline coefficients in fnCreateParam()

# eppspectrum 0.1.2
- Parse Spectrum births in read_specdp_demog_param().
- Add births from Spectrum to specfp object.
- For epp regions, distribute number of births based on 15-44 population size and GFR.
- Implement age 15 entrant prevalence and ART coverage inputs.

# eppspectrum 0.1.1

- Add option to use .ep5 file inputs instead of .DP file when preparing fit from PJNZ file. (use_ep5=FALSE by default.)
- update prepare_spec_fit() to not use popadjust if EPP projection has a single population.


# eppspectrum 0.1

- Update to interface with `epp v0.3` which implements ANC routine testing (ANC-RT) likelihood in EPP 2017.
