
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eppasm

The goal of eppasm is to estimate HIV incidence, prevalence and
mortality using antenatal care data and/or household surveys.

# key populations

This branch can be used for estimation in concentrated epidemics and key
populations (rather than the main branch for urban/rural splits only).
The main differences are i. generalization of PJNZ file extraction code
to work with non-urban/rural split files, ii. incorporation of turnover
into simmod function, iii) adjustments throughout to work with 0
populations (e.g. males for female sex worker population) and alternate
subpopulation names

# remaining issues

While the code works for the example below, there is additional
debugging required to resolve i) currently extremely low incidence (and
therefore prevalence) in key populations, ii) ensure generalization
across more countries, iii) choose more robust base population and entry
rates over time and iv) translate the code to C for speed and efficiency

## examples of current problems using Senegal test file

\[Follow eppasm-basic-key-pops-example.RmD to regenerate results\]

``` r
print(system.file("extdata/example_background_deaths.png", package = "eppasm"))
#> [1] ""
# ![caption.]("inst/extdata/example_ART.png")
# !(system.file("extdata/example_background_deaths.png", package = "eppasm"))
# !(system.file("extdata/example_PLHIV.png", package = "eppasm"))
# !(system.file("extdata/example_HIV_deaths.png", package = "eppasm"))
```
