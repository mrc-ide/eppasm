# Read Spectrum \_pop1.xlsx export file

Spectrum has a debug mode option to output a file recording the entire
state space when a model simulation is computed. This function parse the
\*\_pop1.xlsx file into long-format data.frame.

## Usage

``` r
read_pop1(pop1file, country, years = 2000:2021)
```

## Arguments

- pop1file:

  path to \*\_pop1.xlsx export file

- country:

  string giving the identifying country / region name to be added as a
  column in the output file

- years:

  integer vector of years to extract

## Value

A data.frame consisting of pop1 file in long format

## Details

TODO: tidy up the dependencies in this function: reshape, pryr, readxl.

TODO: Insert instructions for how to enter Spectrum debug mode.

\`years\` must be in the file as sheets, but this is not currently
checked. It will fail hard. Inserting a check for this would make it
fail softer.
