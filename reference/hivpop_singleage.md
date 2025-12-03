# Convert aggregate HIV population to single year

EPP-ASM tracks the CD4 distribution and ART duration of the HIV
population by coarse age groups \`15-16, 17-19, 20-24, ..., 45-49, 50+\`
for computational efficiency.

## Usage

``` r
hivpop_singleage(mod, ss)
```

## Arguments

- mod:

  EPP-ASM model output of class \`spec\`

- ss:

  Model state space inputs

## Details

\`hivpop_singleage\` converts the coarse age group CD4 distribution and
ART coverage to single year of age counts assuming uniform proportions
in each category within coarse age groups.

## Examples

``` r
pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
fp <- prepare_directincid(pjnz)
mod <- simmod(fp)
hivp1 <- hivpop_singleage(mod, fp$ss)
```
