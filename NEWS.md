## epapsm 0.8.1

* Implement Spectrum Adult ART scalar adjustment. This is a user input that 
  allows the input number on ART to be adjusted by a scalar to account for 
  over/under-reporting of treatment numbers.
	
## eppasm 0.8.0

* Implement Spectrum ART allocation.
  
  There has been a longstanding discrepancy betweeen EPP-ASM and Spectrum in ART allocation.
  For ART allocation by 'expected mortality', EPP-ASM allocated according to mortality by CD4
  and age.
  
  Spectrum allocates ART in a two step process: first, ART is allocated by CD4 category based
  on the 'expected mortality' and 'proportional to eligibility' weight. Second, within 
  CD4 categories, ART is allocated by age solely proportional to number in each age 
  group (propotional to eligibility).
  
  This has modest overall difference, but was a source of numerical differences between 
  Spectrum and EPP-ASM.

* Patch ART dropout implementation. Spectrum converts input ART dropout percent to an 
  annual rate using [dropout rate] = -log(1.0 - [input percent]).
  
## eppasm 0.7.7

* Update to use full names for R internal functions e.g. `Rf_allocVector` instead of `allocVector`. Shorthand names are no longer allowed in R v4.5.0. See Nov 10th news https://developer.r-project.org/blosxom.cgi/R-devel

## eppasm 0.7.6

* Update internal data country ISO3 list to contain St. Kitts & Nevis and Dominica

## eppasm 0.7.5

* Qualify all package names and add all required packages into Imports section.

## eppasm 0.7.4

* Implement recovery to next higher CD4 category following ART interruption for those on ART greater than one year.

## eppasm 0.7.3

* Bug fix: account for end-year net migration in the ART population in the first year of ART start.

## eppasm 0.7.2

* Add function `read_pop1()` to parse Spectrum `_pop1.xlsx` export file. 

## eppasm 0.7.1

* Use `vroom::vroom()` to read PJNZ files; much faster for reading `.DP` file.
* Handle Spectrum version numbers saved on Francophone locale devices: e.g. 6,13 instead of 6.13.
* In `read_specdp_demog_param()`, ensure no zero totals when normalising age-specific fertility 
  distribution and net-migration age distribution.
* Fix calculation for target number on treatment during transition from number to percent input.

Package tidying to address R CMD CHECK warnings:

* Add package imports for `anclik` and `binom`.
* Align function signatures for S3 generics.
* Remove `src/Makevars` with GCC specific `-pedantic` compiler flag.

## eppasm 0.7.0

* In `read_specdp_demog_param()`, normalise the age-specific fertility distribution before disaggregating TFR to ASFR. This fixes small discrepancy between Spectrum and EPP-ASM population projection.

* In `read_specdp_demog_param()`, disaggregate the under-5 net migrations to single-year ages proportional to survival probabilities, to match Spectrum.

* Change model projection to calendar year steps instead of mid-year steps, consistent with Spectrum 6.2 updated for WPP 2022 release in December 2022.
  - Net-migrations added at end of projection step, consistent with WPP 2022. No longer (1) adjust net migration for half-period survival, nor (2) adjust 
    net migration to be half in current age group and half in next age group.
  - ART interpolation is for calendar year instead of mid-year to mid-year.
  - Incidence rate 15-49 calculation uses interpolated mid-year susceptible population as denominator.
  - Likelihood calculation use mid-year HIV prevalence or ART coverage based on simple average of end-year values.
  
Code changes are backwards compatible such that code supports simulation of either calendar year or mid-year projection steps. This is controlled by option `fp$projection_period` = either `"midyear"` or `"calendar"`. 

The function `create_spectrum_fixpar(..., projection_period = "calendar")` has new argument which defaults to calendar year input. This must be re-run to recreate the `fp` object to change the choice because calculation of paediatric net migration inputs depends on the choice.


## eppasm 0.6.2

* Distribute age incidence rate ratio using **current** year HIV population by time step, instead of previous year HIV population to match Spectrum calculation. 
  - Calculation of incidence rate by sex uses **previous** year HIV negative population
  - Thanks Rob Glaubius for debugging: https://github.com/mrc-ide/leapfrog/issues/18

## eppasm 0.6.1

* For direct incidence input, create new option `DIRECTINCID_HTS` for intercalating direct incidence new HIV infections in each time step. This was implemented by (1) calculating the number of new infections once per year before the HIV simulation, then (2) adding DT * new infections in each time step. 

Pre-calculating the number of new infections, rather than applying the HIV incidence rate in each time step was preferred to maintain consistency with Spectrum which calculates the number of infections based on previous year population size.

## eppasm 0.6.0

* Add likelihood for survey ART coverage.

## eppasm 0.5.12

* Add condition in `eppasm.cpp` to avoid dividing by 0 when there is zero population eligible
  to initiate ART.

## eppasm 0.5.11

* Patch `get_dp_version()` to search for tag instead of relying on exact location.

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
