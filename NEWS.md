# eppspectrum 0.1

- Update to interface with `epp v0.3` which implements ANC routine testing (ANC-RT) likelihood in EPP 2017.


# eppspectrum 0.1.1

- Add option to use .ep5 file inputs instead of .DP file when preparing fit from PJNZ file. (use_ep5=FALSE by default.)
- update prepare_spec_fit() to not use popadjust if EPP projection has a single population.

# eppspectrum 0.1.3

- EPP model setup occurs in fitmod()
- Allow order of spline model penalty to vary (fp$rtpenord)
- Allow number of knots to vary (fp$numKnots)
- Return spline coefficients in fnCreateParam()

# eppspectrum 0.1.5

- Implemented O-spline model (penalized B-spline) using mgcv to construct
- Added `r0logiotaratio` option to estimate ratio of r(t0) and log of initial seed, reduces correation. Not set as default.

# eppspectrum 0.1.6

- Bug fix: add check in ll_hhage_binom() that ldbinom() does not return NaN.
- Bug fix: in create_spectrum_fixpar(), handle case where ART initiation starts in different years for men and women.


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

# eppasm 0.3.2

- Integrate random effects variance parameter out of likelihood.
- Implement hybrid r-spline and RW model for r(t).
- Implement logistic model for log r(t) (rlogistic), and version with random walk (rlogistic_rw).
- Add optimization option to fitmod().

- Add option to fit without site-level ANC data.

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


# eppasm 0.4.0
- Add annual direct incidence input option.

# eppasm 0.4.1
- Model fitting for case surveillance and vital registration data.
