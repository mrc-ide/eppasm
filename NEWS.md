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

- Specify knot locations for o-spline model

# eppspectrum 0.1.6

- Estimate HIV+:HIV- FRR adjustment.
- Specify knot locations for o-spline model

- Bug fix: add check in ll_hhage_binom() that ldbinom() does not return NaN.
- Bug fix: in create_spectrum_fixpar(), handle case where ART initiation starts in different years for men and women.

# eppspectrum 0.2.0

- Estimate non-HIV mortality trend.
- Add outputs and plot functions for 45q15, 35q15, and age specific mortality.


# eppspectrum 0.2.1

- Implement Wilmoth et al. log-quadratic model.
- Implement incidence rate ratio time x age interaction.
- Trend in h parameter of LogQuad model specified as penalised B-spline with fixed penalty.
- Estimate TIPS coefficients using informative prior distribution.

# eppasm 0.3.0

- Rename package `eppasm`.
- Add IMIS function to eppasm package. In this version of IMIS, at each iteration a new mixture component is constructed either centered on input with greatest weight or based on optimizer. Allows optimizers to be run at arbitrary iterations, but doesn't implment multiple optimizers which might still be usful. Covariance for optimizer mixture component when Hessian is degenerate still needs further work.

- Change from dependency to imports epp package.

# eppasm 0.3.1

- Merge mortality estimation from v0.2.1
