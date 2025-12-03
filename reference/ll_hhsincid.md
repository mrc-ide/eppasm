# Log-likelhood for direct incidence estimate from household survey

Calculate log-likelihood for nationally representative incidence
estimates from a household survey. Currently implements likelihood for a
log-transformed direct incidence estimate and standard error. Needs to
be updated to handle incidence assay outputs.

## Usage

``` r
ll_hhsincid(mod, fp, hhsincid.dat)
```

## Arguments

- mod:

  model output, object of class \`spec\`.

- hhsincid.dat:

  prepared houshold survey incidence estimates (see perp
