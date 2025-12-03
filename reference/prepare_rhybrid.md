# Setup the r-hybrid model

Setup the r-hybrid model

## Usage

``` r
prepare_rhybrid(
  fp,
  tsEpidemicStart = fp$ss$time_epi_start + 0.5,
  rw_start = fp$rw_start,
  rw_trans = fp$rw_trans,
  rw_dk = fp$rw_dk
)
```

## Arguments

- fp:

  model parameters object

- tsEpidemicStart:

  time step at which epidemic is seeded

- rw_start:

  time when random walk starts

- rw_trans:

  number of years to transition from logistic differences to RW
  differences. If NULL, defaults to 5 steps
