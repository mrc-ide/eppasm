# Age-specific ART coverage among pregnant women

Age-specific ART coverage among pregnant women

## Usage

``` r
agepregartcov(
  mod,
  fp,
  aidx = 3:9 * 5 - fp$ss$AGE_START + 1L,
  yidx = 1:fp$ss$PROJ_YEARS,
  agspan = 5,
  expand = FALSE
)
```

## Arguments

- expand:

  whether to expand aidx, yidx, sidx, and agspan
