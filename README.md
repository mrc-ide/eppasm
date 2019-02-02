EPP-ASM
================
2019-02-02

This vignette documents the implementation of age- and risk-group-mixing
(mixing for short) model (Tim’s) into current EPP-ASM code. The code
lives at `mixing` branch. The differences compare to the old version
are:

  - We need a new list named `mx` (perhaps put into `fp` later)
    containing all mixing parameters.
  - Adjustments are made to allow more risk groups (e.g., sex indexing,
    distributing of population into groups)
  - The calculation of the force of infection (FOI), consequently the
    number of new infection in each time step.

The following assumptions and parameters are used for prototyping the
model.

## Force of infection

The force of infection (FOI) is defined as in [Hallett, T. B., Gregson,
S., et al. (2007)](http://doi.org/10.1136/sti.2006.023606). Briefly, the
FOI (![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda
"\\lambda")) of a population of sex
![S=\\{1,2\\}](https://latex.codecogs.com/png.latex?S%3D%5C%7B1%2C2%5C%7D
"S=\\{1,2\\}"), age
![A=1,\\dots](https://latex.codecogs.com/png.latex?A%3D1%2C%5Cdots
"A=1,\\dots"), and risk group ![G
= 1\\dots](https://latex.codecogs.com/png.latex?G%20%3D%201%5Cdots
"G = 1\\dots") is calculated at each time steps as   
![\\lambda\_{SAG} = C\_{SAG}\\sum\_g\\sum\_a C\_{SAGsag} P\_{sag}
(1-K\_{Aa}),](https://latex.codecogs.com/png.latex?%5Clambda_%7BSAG%7D%20%3D%20C_%7BSAG%7D%5Csum_g%5Csum_a%20C_%7BSAGsag%7D%20P_%7Bsag%7D%20%281-K_%7BAa%7D%29%2C
"\\lambda_{SAG} = C_{SAG}\\sum_g\\sum_a C_{SAGsag} P_{sag} (1-K_{Aa}),")  
where ![s,a,g](https://latex.codecogs.com/png.latex?s%2Ca%2Cg "s,a,g")
denote the sex, age, and risk group of the opposite sex, respectively.
![C\_{SAG}](https://latex.codecogs.com/png.latex?C_%7BSAG%7D "C_{SAG}")
is the total number of contacts,
![C\_{SAGsag}](https://latex.codecogs.com/png.latex?C_%7BSAGsag%7D
"C_{SAGsag}") is the mixing pattern with other age and group,
![P\_{sag}](https://latex.codecogs.com/png.latex?P_%7Bsag%7D "P_{sag}")
is the tranmission probability, and
![K\_{Aa}](https://latex.codecogs.com/png.latex?K_%7BAa%7D "K_{Aa}") is
the probability of correct and consistent condom usage given the pair.
All the model’s terms are time-dependent but
![K\_{Aa}](https://latex.codecogs.com/png.latex?K_%7BAa%7D "K_{Aa}")
might not, depending on the availability of data.

### Risk groups mixing

![C\_{sag}](https://latex.codecogs.com/png.latex?C_%7Bsag%7D "C_{sag}")
is calculated as   
![C\_{SAG} = \\dfrac{M}{\\tau^{\\sum m\\gamma\_m}}\\tau^m, \\qquad
m=G-1,
](https://latex.codecogs.com/png.latex?C_%7BSAG%7D%20%3D%20%5Cdfrac%7BM%7D%7B%5Ctau%5E%7B%5Csum%20m%5Cgamma_m%7D%7D%5Ctau%5Em%2C%20%5Cqquad%20m%3DG-1%2C%20
"C_{SAG} = \\dfrac{M}{\\tau^{\\sum m\\gamma_m}}\\tau^m, \\qquad m=G-1, ")  
where ![M](https://latex.codecogs.com/png.latex?M "M") is the geometric
mean of the total number of contacts,
![\\tau](https://latex.codecogs.com/png.latex?%5Ctau "\\tau") is common
ratio between risk groups and
![\\gamma\_m](https://latex.codecogs.com/png.latex?%5Cgamma_m
"\\gamma_m") is the proportion of population in the group
![G](https://latex.codecogs.com/png.latex?G "G").

The proportion in low and high risk group are assumed similar for both
male and female. This is used to distribute:

  - the initial population
  - the entrants, healthy and H+
  - immigrants

> TODO: allow different for male and female

``` r
mx <- list()
mx$gamma <- c(.9, .1) # order Low -> High Risk, sum to 1
```

If `mx$gamma` input is of length one, then it is set to 1 and one risk
group model is run.

Higher risk group have tripple the frequency of contacts of the
immediate lower risk group.

``` r
mx$tau <- 3 # common ratio between risk groups
```

where the number of contacts is given with an age- and sex-specific
geometric mean number of contacts.

``` r
mx$M <- matrix(c(c(1.6, 1.6, 1.6, 1.7, 1.6, 1),
                 c(1.2, 1.2, 1.1, 1.5, 1.1, 1),
                 c(15, 20, 25, 30, 40, 50),
                 c(19, 24, 29, 39, 49, 80)), 6, 4, dimnames=list(c("1529", "2024", "2529", "3039", "4049", "50+"), c("M", "F", "lo", "up")))
# expanding the age group provided in M
mx$M. <- function(x) apply(x[, 1:2], 2, rep, times = apply(x[, 3:4], 1, diff) + 1)
mx$Mx <- mx$M.(mx$M)
```

|      |   M |   F | lo | up |
| ---- | --: | --: | -: | -: |
| 1529 | 1.6 | 1.2 | 15 | 19 |
| 2024 | 1.6 | 1.2 | 20 | 24 |
| 2529 | 1.6 | 1.1 | 25 | 29 |
| 3039 | 1.7 | 1.5 | 30 | 39 |
| 4049 | 1.6 | 1.1 | 40 | 49 |
| 50+  | 1.0 | 1.0 | 50 | 80 |

Table: Risk group mixing

![C\_{SAGsag}](https://latex.codecogs.com/png.latex?C_%7BSAGsag%7D
"C_{SAGsag}") is calculated as   
![C\_{SAGsag} = \\left\[ (1-\\epsilon)\\delta\_{Gg} +
\\epsilon\\left(\\dfrac{C\_{sag}N\_{sag}}{\\sum\_g N\_{sag}
C\_{sag}}\\right)\\right\]\\Delta\_{SAsa},](https://latex.codecogs.com/png.latex?C_%7BSAGsag%7D%20%3D%20%5Cleft%5B%20%281-%5Cepsilon%29%5Cdelta_%7BGg%7D%20%2B%20%5Cepsilon%5Cleft%28%5Cdfrac%7BC_%7Bsag%7DN_%7Bsag%7D%7D%7B%5Csum_g%20N_%7Bsag%7D%20C_%7Bsag%7D%7D%5Cright%29%5Cright%5D%5CDelta_%7BSAsa%7D%2C
"C_{SAGsag} = \\left[ (1-\\epsilon)\\delta_{Gg} + \\epsilon\\left(\\dfrac{C_{sag}N_{sag}}{\\sum_g N_{sag} C_{sag}}\\right)\\right]\\Delta_{SAsa},")  
where ![\\epsilon \\in
\[0,1\]](https://latex.codecogs.com/png.latex?%5Cepsilon%20%5Cin%20%5B0%2C1%5D
"\\epsilon \\in [0,1]") defines whether contacts are assortative or
random with
![\\delta\_{Gg}](https://latex.codecogs.com/png.latex?%5Cdelta_%7BGg%7D
"\\delta_{Gg}") is the Koronecker delta.
![N\_{sag}](https://latex.codecogs.com/png.latex?N_%7Bsag%7D "N_{sag}")
is the size of the corresponding population.
![\\Delta\_{SAsa}](https://latex.codecogs.com/png.latex?%5CDelta_%7BSAsa%7D
"\\Delta_{SAsa}") is the age-mixing pattern of the pair, defining the
probability that a partnership is formed between an aged
![A](https://latex.codecogs.com/png.latex?A "A") subject and aged
![a](https://latex.codecogs.com/png.latex?a "a") subject.

``` r
mx$epsilon <- .2
```

If there is only one group,
![\\epsilon](https://latex.codecogs.com/png.latex?%5Cepsilon
"\\epsilon") will be set to zero by checking the length of
![\\gamma](https://latex.codecogs.com/png.latex?%5Cgamma "\\gamma")
provided.

The following rates are the same among risk groups, including `Sx,
cd4_mort, cd4_initdist, cd4_prog, paedsurv_cd4dist, art_mort,
paedsurv_artcd4dist`.

![C\_{SAGsag}](https://latex.codecogs.com/png.latex?C_%7BSAGsag%7D
"C_{SAGsag}") is adjusted to ensure the balance of the number of
partnership calculated from each group, i.e.,

  
![N\_{SAG} C\_{SAGsag} = N\_{sag}
C\_{sagSAG}.](https://latex.codecogs.com/png.latex?N_%7BSAG%7D%20C_%7BSAGsag%7D%20%3D%20N_%7Bsag%7D%20C_%7BsagSAG%7D.
"N_{SAG} C_{SAGsag} = N_{sag} C_{sagSAG}.")  

In particular, we calculate   
![W\_{SAGsag} = \\dfrac{N\_{SAG} C\_{SAGsag}}{N\_{sag}
C\_{sagSAG}}](https://latex.codecogs.com/png.latex?W_%7BSAGsag%7D%20%3D%20%5Cdfrac%7BN_%7BSAG%7D%20C_%7BSAGsag%7D%7D%7BN_%7Bsag%7D%20C_%7BsagSAG%7D%7D
"W_{SAGsag} = \\dfrac{N_{SAG} C_{SAGsag}}{N_{sag} C_{sagSAG}}")  
and adjust

  
![
C\_{SAGsag} = C\_{SAGsag} \\sqrt{W\_{SAGsag}}
](https://latex.codecogs.com/png.latex?%0A%20%20C_%7BSAGsag%7D%20%3D%20C_%7BSAGsag%7D%20%5Csqrt%7BW_%7BSAGsag%7D%7D%0A
"
  C_{SAGsag} = C_{SAGsag} \\sqrt{W_{SAGsag}}
")  

  
![
C\_{sagSAG} &= C\_{sagSAG} / \\sqrt{W\_{SAGsag}}
](https://latex.codecogs.com/png.latex?%0A%20%20C_%7BsagSAG%7D%20%26%3D%20C_%7BsagSAG%7D%20%2F%20%5Csqrt%7BW_%7BSAGsag%7D%7D%0A
"
  C_{sagSAG} &= C_{sagSAG} / \\sqrt{W_{SAGsag}}
")  

In total, we calculate ![2\\times2\\times
G^2](https://latex.codecogs.com/png.latex?2%5Ctimes2%5Ctimes%20G%5E2
"2\\times2\\times G^2") matrices for this weighting process.

The risk group membership is fixed, i.e. no moving from one risk group
to another.

### Age mixing

![\\Delta\_{Aa}](https://latex.codecogs.com/png.latex?%5CDelta_%7BAa%7D
"\\Delta_{Aa}") can be defined similarly to the group mixing term with
another
![\\varepsilon](https://latex.codecogs.com/png.latex?%5Cvarepsilon
"\\varepsilon") defining the age-mixing pattern. Here
![\\Delta\_{Aa}](https://latex.codecogs.com/png.latex?%5CDelta_%7BAa%7D
"\\Delta_{Aa}") is defined by fitting a log-logistic model
(![f](https://latex.codecogs.com/png.latex?f "f")) to the sexual partner
age difference data. For female aged
![A](https://latex.codecogs.com/png.latex?A "A") and male aged
![a](https://latex.codecogs.com/png.latex?a "a")

  
![\\Delta\_{Aa} = 
\\begin{cases}
f(a -A),& \\qquad a \\geq A\\\\
0,& \\text{otherwise},
\\end{cases}
](https://latex.codecogs.com/png.latex?%5CDelta_%7BAa%7D%20%3D%20%0A%5Cbegin%7Bcases%7D%0A%20%20f%28a%20-A%29%2C%26%20%5Cqquad%20a%20%5Cgeq%20A%5C%5C%0A%20%200%2C%26%20%5Ctext%7Botherwise%7D%2C%0A%5Cend%7Bcases%7D%0A
"\\Delta_{Aa} = 
\\begin{cases}
  f(a -A),& \\qquad a \\geq A\\\\
  0,& \\text{otherwise},
\\end{cases}
")  
and for male ![\\Delta\_{Aa} =
\\Delta\_{aA}](https://latex.codecogs.com/png.latex?%5CDelta_%7BAa%7D%20%3D%20%5CDelta_%7BaA%7D
"\\Delta_{Aa} = \\Delta_{aA}"). This is adjusted such that
![\\sum\_{a=a\_L}^{a=a\_U}\\Delta\_{Aa}
= 1,](https://latex.codecogs.com/png.latex?%5Csum_%7Ba%3Da_L%7D%5E%7Ba%3Da_U%7D%5CDelta_%7BAa%7D%20%3D%201%2C
"\\sum_{a=a_L}^{a=a_U}\\Delta_{Aa} = 1,") where ![a\_L,
a\_U](https://latex.codecogs.com/png.latex?a_L%2C%20a_U "a_L, a_U")
denote the age at first sex and age at last sex, respectively.

Age- and sex-specific mixing matrix
(![\\Delta\_{Aa}](https://latex.codecogs.com/png.latex?%5CDelta_%7BAa%7D
"\\Delta_{Aa}")) follows a three parameters log-logistic distribution
(![\\kappa, \\rho,
r](https://latex.codecogs.com/png.latex?%5Ckappa%2C%20%5Crho%2C%20r
"\\kappa, \\rho, r")) in which the scale parameter `r` is shifted by the
age difference between male and female.

``` r
mx$kappa <- 2 # log-logistic shape
mx$rho <- 0.25 # log-logistic scale
r <- 0.15 # log-logistic support
```

This distribution seems too strict. Considering options that allows
female with younger male and higher probability for the same age.

![Fig: Age mixing](README_files/figure-gfm/unnamed-chunk-7-1.png)

It might be easier to work with the normal age-mixing formula with only
one parameter, but need data to fit.

### Condom usage

The probability of correct and consistent use of among age-groups and
sexes are given as inputs. Below is a condom use data for only ±25 years
old. This is age-dependent not risk group dependent which should be
modified depend on data availability. It is assumed that condom efficacy
in preventing HIV transmission is 90%.

``` r
# prob of correct and consistent condom use by <25 and 25+
mx$chi <- matrix(c(.32, .1, .39, .07), 2, 2,
                 dimnames=list(c("M<25", "M25+"), c("F<25", "F25+")))
mx$cEf <- .9
```

|       | F\<25 | F25+ |
| ----- | ----: | ---: |
| M\<25 |  0.32 | 0.39 |
| M25+  |  0.10 | 0.07 |

Table: Condom usage

### Sexual debut

Population under the age of first sex and above the age of last sex
would not have any sexual contacts in this model. For now the code
covers only the case that both sexes have the same sexual debut period.

``` r
mx$A1 = 16 # age at first sex
mx$A2 = 55 # age at last sex
```

Setting this effectively set FOI of population outside the boundary zero
by making ![\\Delta\_{AA}
= 0](https://latex.codecogs.com/png.latex?%5CDelta_%7BAA%7D%20%3D%200
"\\Delta_{AA} = 0") for those not having sexual acts.

> TODO: different for risk groups

### Tranmission rate

Transmission rate is assumed different for each CD4+ stages (![d
= 1,\\dots,7](https://latex.codecogs.com/png.latex?d%20%3D%201%2C%5Cdots%2C7
"d = 1,\\dots,7")) and on ART or not (![k
= 1, 0](https://latex.codecogs.com/png.latex?k%20%3D%201%2C%200
"k = 1, 0")).

  
![
P\_{sag} = r\_t \\cdot
\\dfrac{\\sum\_{d} w\_{d} \[H\_{sag}^{d,k=0} - H\_{sag}^{d,k=1} (1 -
\\phi)\]}{N\_{sag}} 
\= r\_t \\cdot \\dfrac{H^{\*}\_{sag}}{N\_{sag}}
](https://latex.codecogs.com/png.latex?%0AP_%7Bsag%7D%20%3D%20r_t%20%5Ccdot%0A%5Cdfrac%7B%5Csum_%7Bd%7D%20w_%7Bd%7D%20%5BH_%7Bsag%7D%5E%7Bd%2Ck%3D0%7D%20-%20H_%7Bsag%7D%5E%7Bd%2Ck%3D1%7D%20%281%20-%20%5Cphi%29%5D%7D%7BN_%7Bsag%7D%7D%20%0A%3D%20r_t%20%5Ccdot%20%5Cdfrac%7BH%5E%7B%2A%7D_%7Bsag%7D%7D%7BN_%7Bsag%7D%7D%0A
"
P_{sag} = r_t \\cdot
\\dfrac{\\sum_{d} w_{d} [H_{sag}^{d,k=0} - H_{sag}^{d,k=1} (1 - \\phi)]}{N_{sag}} 
= r_t \\cdot \\dfrac{H^{*}_{sag}}{N_{sag}}
")  

where ![H](https://latex.codecogs.com/png.latex?H "H") denotes the HIV
positive population, ![w\_d](https://latex.codecogs.com/png.latex?w_d
"w_d") denotes the relative infectiousness of individuals in the disease
stage ![d](https://latex.codecogs.com/png.latex?d "d"),
![\\phi](https://latex.codecogs.com/png.latex?%5Cphi "\\phi") is the
relative reduction in infectiousness of the individual taking
cART.

``` r
mx$wD <- cbind(c( 500, 350, 250, 200, 100,  50,  0), # lower bound CD4+ count
               c(1000, 500, 350, 250, 200, 100, 50), # upper
               c(   1,  .1,  .2,  .2,  .7,  .7, .7)) # relative infectiousness
mx$phi <- .3 # ~fp$relinfectART
```

But as the HIV population is grouped,
![H^{\*}\_{sag}](https://latex.codecogs.com/png.latex?H%5E%7B%2A%7D_%7Bsag%7D
"H^{*}_{sag}") is calculated by age-group and distributed as   
![H^{\*}\_{sag} = \\dfrac{N\_{sag}}{N\_{s\\alpha g}} H^{\*}\_{s\\alpha
g}
](https://latex.codecogs.com/png.latex?H%5E%7B%2A%7D_%7Bsag%7D%20%3D%20%5Cdfrac%7BN_%7Bsag%7D%7D%7BN_%7Bs%5Calpha%20g%7D%7D%20H%5E%7B%2A%7D_%7Bs%5Calpha%20g%7D%20
"H^{*}_{sag} = \\dfrac{N_{sag}}{N_{s\\alpha g}} H^{*}_{s\\alpha g} ")  

where ![\\alpha](https://latex.codecogs.com/png.latex?%5Calpha
"\\alpha") the age-group contains age
![a](https://latex.codecogs.com/png.latex?a "a").

## Model simulation

Prepare the fix parameters as the vignette example

### Without mixing

Simulate the old model, to check if the tidy up process and the added
code do break anything.

``` r
devtools::load_all()
#> Loading eppasm
modOC <- simmod(fp_par, "C") # old C
#> 
#> Attaching package: 'magrittr'
#> The following objects are masked from 'package:testthat':
#> 
#>     equals, is_less_than, not
modOR <- simmod(fp_par, "R") # old R
prev(modOC)
#>  [1] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
#>  [6] 0.0000000000 0.0002115617 0.0003353572 0.0005295253 0.0008332484
#> [11] 0.0013069493 0.0020358322 0.0031619043 0.0048919274 0.0075278060
#> [16] 0.0114952631 0.0173160779 0.0256766875 0.0372965514 0.0527891555
#> [21] 0.0724249143 0.0959981021 0.1223878907 0.1498574784 0.1763411128
#> [26] 0.1999167706 0.2192912261 0.2337275925 0.2432604836 0.2483679227
#> [31] 0.2497462402 0.2481788050 0.2442821044 0.2390433031 0.2336087314
#> [36] 0.2289475455 0.2248869398 0.2212650129 0.2179111213 0.2148851567
#> [41] 0.2120610141 0.2092697998 0.2057857477 0.2021365230 0.1981385444
#> [46] 0.1934484450 0.1879411839 0.1818498064 0.1757016454 0.1701286410
#> [51] 0.1646202991 0.1589677009
prev(modOR)
#>  [1] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
#>  [6] 0.0000000000 0.0002115617 0.0003353572 0.0005295253 0.0008332483
#> [11] 0.0013069491 0.0020358315 0.0031619027 0.0048919239 0.0075277994
#> [16] 0.0114952515 0.0173160589 0.0256766580 0.0372965076 0.0527890932
#> [21] 0.0724248298 0.0959979926 0.1223877553 0.1498573185 0.1763409318
#> [26] 0.1999165733 0.2192910176 0.2337273776 0.2432602664 0.2483677063
#> [31] 0.2497460271 0.2481785968 0.2442819025 0.2390431089 0.2336085461
#> [36] 0.2289473692 0.2248867720 0.2212648522 0.2179109669 0.2148850069
#> [41] 0.2120608674 0.2092696564 0.2057856066 0.2021363852 0.1981384097
#> [46] 0.1934483116 0.1879410523 0.1818496774 0.1757015190 0.1701285554
#> [51] 0.1646202279 0.1589676418
```

The later years are slightly differrence after the third digit. This
results from correcting the negative target population inputs and remove
one line on `grad` update (?). These are important for the mixing model
as the sign of population trouble FOI calculation.

### Mixing model

Model with one and two risk groups per sex.

``` r
devtools::load_all() # 2b remove
#> Loading eppasm
mx$gamma <- c(.9,.1)
modN2 <- simmod(fp_par, "R", isMixing = TRUE, mx)
mx$gamma <- 1
modN1 <- simmod(fp_par, "R", isMixing = TRUE, mx)
# profvis::profvis(simmod(fp_par, "R", isMixing = TRUE, mx))
# Vectorize outer make thing worse
```

FOI

``` r
plot(modN2)
```

![Figure: FOI two risk
groups](README_files/figure-gfm/unnamed-chunk-15-1.png)

``` r
plot(modN1)
```

![Figure: FOI two risk
groups](README_files/figure-gfm/unnamed-chunk-16-1.png)

Prevalence

![Figure: FOI two risk
groups](README_files/figure-gfm/unnamed-chunk-17-1.png)

> TODO: add relavtive transmission rate for pair of risk group. This
> should go to transmission rate section.

``` r
# Template data (Thembisa)
mx$betaPair <- array(1, c(2,2,2))
# Susceptible male
mx$betaPair[,,1] <- matrix(c(1, 1.25, 1.08, 1.23), 2, 2,
                      dimnames=list(c("M-Low", "M-High"), c("F-Low", "F-High")))
# Susceptible female
mx$betaPair[,,2] <- matrix(c(1, 1.14, 1.09, 1.20), 2, 2,
                      dimnames=list(c("M-Low", "M-High"), c("F-Low", "F-High")))
```

|      |      |
| ---: | ---: |
| 1.00 | 1.08 |
| 1.25 | 1.23 |

Table: Relative transmission rate for susceptible male

|      |      |
| ---: | ---: |
| 1.00 | 1.09 |
| 1.14 | 1.20 |

Table: Relative transmission rate for susceptible female

> TODO: which data, fitting which parameteres? TODO: Projection?