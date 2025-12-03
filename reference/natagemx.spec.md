# Non-HIV age-specific mortality

Calculate all-cause mortality rate by single year of age and sex from a
`spec` object.

## Usage

``` r
# S3 method for class 'spec'
natagemx(mod)
```

## Arguments

- mod:

  output of simmod of class
  [`spec`](https://rdrr.io/r/stats/spectrum.html).

## Value

3-dimensional array of mortality by age, sex, and year.

## Details

Mortality in year Y is calculated as the number of non-HIV deaths
occurring from the mid-year of year Y-1 to mid-year Y, divided by the
population size at the mid-year of year Y-1. !!! NOTE: This might be
different from the calculation in Spectrum. Should confirm this with
John Stover.
