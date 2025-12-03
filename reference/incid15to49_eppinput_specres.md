# Incidence rate age 15-49 using previous-year susceptible population

Return HIV incidence calculated as number of infections among age 15-49
population divided by susceptible population at the start of the
projection year. This was the standard HIV incidence calculation
reported by Spectrum up to version \<=6.19. From Spectrum version
\>=6.2, the incidence rate is reported divided by the projection period
population.

## Usage

``` r
incid15to49_eppinput_specres(x)
```

## Arguments

- x:

  \`specres\` object created by \`read_hivproj_output()\`.
