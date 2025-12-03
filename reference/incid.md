# Incidence rate among adults age 15-49 years

Incidence rate among adults age 15-49 years

## Usage

``` r
incid(mod, ...)
```

## Arguments

- mod:

  model output

## Details

This returns incidence rate calculated as the number of infections
during the projection period divded by the number susceptible at the
mid-point of the projection period. This is the default incidence
calculation for Spectrum version \>=6.2 For Spectrum versions \<=6.19
incidence was calculated as number of infections divided by susceptible
population at the start of the projection year.
