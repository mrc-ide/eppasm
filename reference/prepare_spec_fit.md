# Prepare fit by EPP regions

Prepare fit by EPP regions

## Usage

``` r
prepare_spec_fit(
  pjnz,
  proj.end = 2016.5,
  popadjust = NULL,
  popupdate = TRUE,
  use_ep5 = FALSE
)
```

## Arguments

- pjnz:

  file path to Spectrum PJNZ file.

- proj.end:

  end year for projection.

- popupdate:

  logical should target population be updated to match age-specific
  population size from DP file and %Urban from EPP XML.
