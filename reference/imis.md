# Incremental Mixture Importance Sampling (IMIS)

Implements IMIS algorithm with optional optimization step (Raftery and
Bao 2010).

## Usage

``` r
imis(
  B0,
  B,
  B_re,
  number_k,
  opt_k = NULL,
  fp,
  likdat,
  prior = eppasm::prior,
  likelihood = eppasm::likelihood,
  sample_prior = eppasm::sample.prior,
  dsamp = eppasm::dsamp,
  save_all = FALSE
)
```

## Arguments

- B0:

  number of initial samples to draw

- B:

  number of samples at each IMIS iteration

- B_re:

  number of resamples

- number_k:

  maximum number of iterations

- opt_k:

  vector of iterations at which to use optimization step to identify new
  mixture component

- fp:

  fixed model parameters

- likdat:

  likeihood data

- prior:

  function to calculate prior density for matrix of parameter inputs

- likelihood:

  function to calculate likelihood for matrix of parameter inputs

- sample_prior:

  function to draw an initial sample of parameter inputs

- dsamp:

  function to calculate density for initial sampling distribution (may
  be equal to prior)

- save_all:

  logical whether to save all sampled parameters

## Value

list with items resample, stat, and center
