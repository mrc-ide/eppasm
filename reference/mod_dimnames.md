# Add dimnames to spec model output

\`mod_dimnames\` assigns dimnames to spec model outputs.

## Usage

``` r
mod_dimnames(mod, ss)
```

## Arguments

- mod:

  EPP-ASM model output of class \`spec\`

- ss:

  Model state space inputs

## Value

EPP-ASM model output with labelled dimensions

## Examples

``` r
pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
fp <- prepare_directincid(pjnz)
mod <- simmod(fp)
```
