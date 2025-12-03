# Create simulation inputs fixed parameters

Create simulation inputs fixed parameters

## Usage

``` r
create_spectrum_fixpar(
  projp,
  demp,
  hiv_steps_per_year = 10L,
  proj_start = projp$yr_start,
  proj_end = projp$yr_end,
  AGE_START = 15L,
  relinfectART = projp$relinfectART,
  time_epi_start = projp$t0,
  popadjust = FALSE,
  targetpop = demp$basepop,
  artelig200adj = TRUE,
  who34percelig = 0,
  frr_art6mos = projp$frr_art6mos,
  frr_art1yr = projp$frr_art6mos,
  projection_period = NULL,
  art_dropout_recover_cd4 = NULL
)
```

## Details

If argument \`projection_period = NULL\`, R determines the projection
period based on the Spectrum version number. For version \<= 6.19,
projection period is \`"midyear"\`, and for version \>= 6.20, projection
period is \`"calendar"\`.
