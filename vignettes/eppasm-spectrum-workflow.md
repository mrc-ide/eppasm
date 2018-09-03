---
title: "EPP-ASM to Spectrum inputs workflow"
author: "Jeff Eaton, Rob Glaubius, Tim Brown, John Stover"
date: "2018-06-29"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



This vignette develops the workflow for fitting the EPP-ASM model and generating outputs to provide as Spectrum inputs.


```r
## eppasm version >= 0.5.2
## devtool::install_github("mrc-ide/eppasm@dev")
## library(eppasm)
devtools::load_all()
#> Loading eppasm
```

## Preparing EPP-ASM inputs


```r
pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
bw <- prepare_spec_fit(pjnz, proj.end=2021.5)
```

## Fitting the EPP-ASM model


```r
bwfit <- list()

bwfit$Urban <- fitmod(bw$Urban,
                      eppmod = "rlogistic_rw", rw_start = 2005, fitincrr = "linincrr",
                      B0=1e4, B=1e3, opt_iter = 1:2*5, number_k=50)
#> [1] "Stage   MargLike   UniquePoint   MaxWeight      ESS   IterTime"
#> [1] "    1     56.090          1.00        1.00     1.00      59.23"
#> [1] "    2     38.518          8.22        0.65     2.07       5.93"
#> [1] "    3     37.434          7.95        0.70     1.89       5.86"
#> [1] "    4     37.353          6.00        0.96     1.08       5.90"
#> [1] "    5     35.848          9.25        0.70     1.93       5.87"
#> [1] "maximum log posterior= 124.13 , time used= 1.71 minutes, convergence= 0"
#> [1] "    6     71.905        302.67        0.26    13.28     110.39"
#> [1] "    7     71.413        325.78        0.05    82.47       6.00"
#> [1] "    8     71.275        326.56        0.05    88.92       5.93"
#> [1] "    9     71.166        324.75        0.05    93.19       5.95"
#> [1] "   10     70.994        325.95        0.04   110.76       6.03"
#> [1] "maximum log posterior= 124.11 , time used= 1.69 minutes, convergence= 0"
#> [1] "   11     71.651        531.77        0.08    56.08     109.24"
#> [1] "   12     71.463        555.89        0.07    83.61       6.09"
#> [1] "   13     71.296        572.71        0.05   120.01       6.04"
#> [1] "   14     71.214        576.51        0.04   160.18       6.06"
#> [1] "   15     71.140        577.42        0.02   203.36       5.79"
#> [1] "   16     71.081        574.65        0.02   211.19       6.11"
#> [1] "   17     71.014        570.98        0.02   216.03       6.11"
#> [1] "   18     70.943        572.84        0.02   223.97       6.06"
#> [1] "   19     70.894        567.72        0.02   228.37       5.95"
#> [1] "   20     70.833        570.29        0.02   243.15       6.06"
#> [1] "   21     70.799        564.22        0.02   246.00       6.18"
#> [1] "   22     70.757        560.31        0.01   246.15       6.05"
#> [1] "   23     70.717        556.28        0.02   244.33       5.99"
#> [1] "   24     70.654        550.92        0.02   239.63       6.27"
#> [1] "   25     70.593        547.02        0.01   235.30       6.07"
#> [1] "   26     70.512        547.61        0.01   240.42       6.42"
#> [1] "   27     70.493        546.09        0.01   244.16       5.99"
#> [1] "   28     70.467        541.34        0.01   242.46       6.15"
#> [1] "   29     70.428        535.39        0.01   237.30       6.30"
#> [1] "   30     70.393        529.25        0.01   232.57       6.19"
#> [1] "   31     70.360        523.53        0.01   231.18       6.25"
#> [1] "   32     70.340        520.38        0.01   231.47       6.08"
#> [1] "   33     70.308        514.36        0.01   229.28       6.23"
#> [1] "   34     70.273        508.00        0.01   224.23       6.25"
#> [1] "   35     70.227        504.48        0.01   219.31       6.42"
#> [1] "   36     70.208        501.11        0.01   218.85       6.25"
#> [1] "   37     70.178        495.95        0.01   215.29       6.26"
#> [1] "   38     70.138        491.24        0.01   210.25       6.27"
#> [1] "   39     70.117        488.99        0.01   210.26       6.13"
#> [1] "   40     70.094        483.26        0.01   208.40       6.23"
#> [1] "   41     70.063        479.67        0.01   205.26       6.62"
#> [1] "   42     70.041        477.05        0.01   204.27       6.11"
#> [1] "   43     70.024        474.45        0.01   203.98       6.23"
#> [1] "   44     69.973        474.79        0.01   200.90       6.74"
#> [1] "   45     69.942        469.91        0.01   200.19       6.35"
#> [1] "   46     69.918        469.23        0.01   199.48       6.36"
#> [1] "   47     69.893        466.87        0.01   198.65       6.32"
#> [1] "   48     69.838        462.09        0.01   194.19       6.70"
#> [1] "   49     69.805        458.24        0.01   190.87       6.39"
#> [1] "   50     69.752        456.09        0.01   190.34       6.60"
bwfit$Rural <- fitmod(bw$Rural,
                      eppmod = "rlogistic_rw", rw_start = 2005, fitincrr = "linincrr",
                      B0=1e4, B=1e3, opt_iter = 1:2*5, number_k=50)
#> [1] "Stage   MargLike   UniquePoint   MaxWeight      ESS   IterTime"
#> [1] "    1     23.080          3.10        0.99     1.02      56.91"
#> [1] "    2     22.136          3.04        0.99     1.02       5.99"
#> [1] "    3     19.019          6.33        0.83     1.39       6.15"
#> [1] "    4     17.716         10.20        0.57     2.44       6.16"
#> [1] "    5     17.141         11.55        0.43     3.20       5.79"
#> [1] "maximum log posterior= 100.5 , time used= 1.51 minutes, convergence= 0"
#> [1] "    6     43.356         47.05        0.49     2.49      98.08"
#> [1] "    7     40.851        122.99        0.47     4.23       6.48"
#> [1] "    8     39.351        186.56        0.16    15.83       6.82"
#> [1] "    9     38.647        218.20        0.28    11.31       7.07"
#> [1] "   10     37.965        251.64        0.07    41.55       7.16"
#> [1] "maximum log posterior= 101.25 , time used= 1.48 minutes, convergence= 0"
#> [1] "   11     43.892        108.44        0.34     5.58      96.29"
#> [1] "   12     43.454        120.47        0.31     6.39       7.06"
#> [1] "   13     42.985        140.68        0.31     6.55       7.44"
#> [1] "   14     42.315        184.73        0.26     8.94       7.16"
#> [1] "   15     41.912        212.57        0.27    10.78       7.38"
#> [1] "   16     41.581        236.57        0.09    27.37       7.17"
#> [1] "   17     41.408        247.12        0.09    26.92       7.19"
#> [1] "   18     41.282        250.57        0.10    26.62       7.13"
#> [1] "   19     41.089        264.79        0.09    25.87       7.15"
#> [1] "   20     40.991        271.29        0.09    26.72       7.36"
#> [1] "   21     40.870        277.39        0.09    26.77       7.22"
#> [1] "   22     40.687        287.06        0.10    24.71       7.37"
#> [1] "   23     40.553        291.61        0.10    25.92       7.32"
#> [1] "   24     40.418        301.62        0.10    27.26       7.46"
#> [1] "   25     40.211        322.72        0.09    29.82       7.38"
#> [1] "   26     40.070        336.68        0.09    29.85       7.34"
#> [1] "   27     39.940        342.83        0.09    31.55       7.40"
#> [1] "   28     39.824        350.81        0.08    35.04       7.44"
#> [1] "   29     39.724        356.13        0.09    37.95       7.43"
#> [1] "   30     39.623        362.66        0.09    43.48       7.51"
#> [1] "   31     39.477        376.68        0.05    56.48       7.54"
#> [1] "   32     39.406        377.57        0.06    58.94       7.48"
#> [1] "   33     39.333        376.90        0.06    63.00       7.41"
#> [1] "   34     39.266        376.16        0.05    68.95       7.72"
#> [1] "   35     39.191        378.99        0.04    70.92       7.12"
#> [1] "   36     39.125        382.70        0.04    73.12       7.45"
#> [1] "   37     39.064        383.49        0.03    75.93       7.50"
#> [1] "   38     38.959        393.03        0.03    76.28       7.67"
#> [1] "   39     38.880        398.60        0.04    76.90       7.63"
#> [1] "   40     38.832        400.90        0.04    78.11       7.52"
#> [1] "   41     38.775        404.36        0.04    79.70       7.44"
#> [1] "   42     38.722        404.08        0.04    80.90       8.06"
#> [1] "   43     38.669        401.59        0.04    82.96       7.35"
#> [1] "   44     38.627        404.08        0.03    86.56       7.41"
#> [1] "   45     38.550        405.05        0.03    86.04       7.35"
#> [1] "   46     38.515        406.64        0.03    88.48       7.50"
#> [1] "   47     38.476        406.10        0.03    90.45       7.58"
#> [1] "   48     38.399        409.58        0.03    88.94       7.54"
#> [1] "   49     38.362        409.76        0.03    90.78       7.84"
#> [1] "   50     38.316        410.99        0.03    92.21       7.50"
```

When fitting, the random-walk based models only simulate through the end of the
data period. The `extend_projection()` function extends the random walk for $r(t)$
through the end of the projection period.


```r
bwfit <- lapply(bwfit, extend_projection, proj_years = 52)
```

## Pooling EPP subpopulation results

The function `aggr_specfit()` 
This involves simulating the model for all resamples in each subregion and summing the following `pop`, `hivpop`, and `artpop` arrays for each of the 3000 resamples to generate 3000 national outputs.


```r
bwaggr <- aggr_specfit(bwfit)
```

## Generating outputs 

Now generate outputs for the following by age, sex, and year:

* Total population
* HIV positive population
* New infections
* Non-HIV deaths
* HIV deaths

This is all a bit manual right now. In due course, will write better functions to assist with this.



```r
library(magrittr)
library(data.table)

dimnm <- list(age = 15:80, sex = c("male", "female"), year = 1970:2021)

totpop <- lapply(bwaggr, apply, c(1, 2, 4), sum) %>%
  abind::abind(., rev.along=0, use.dnns=TRUE,
               new.names = c(dimnm, list(sampleid = seq_along(.)))) %>%
  melt %>%
  data.table(., outcome = "totpop")

hivpop <- lapply(bwaggr, function(x) x[ , , 2, ]) %>%
  abind::abind(., rev.along=0, use.dnns=TRUE,
               new.names = c(dimnm, list(sampleid = seq_along(.)))) %>%
  melt %>%
  data.table(., outcome = "hivpop")

infections <- lapply(bwaggr, attr, "infections") %>%
  abind::abind(., rev.along=0, use.dnns=TRUE,
               new.names = c(dimnm, list(sampleid = seq_along(.)))) %>%
  melt %>%
  data.table(., outcome = "infections")

natdeaths <- lapply(bwaggr, attr, "natdeaths") %>%
  abind::abind(., rev.along=0, use.dnns=TRUE,
               new.names = c(dimnm, list(sampleid = seq_along(.)))) %>%
  melt %>%
  data.table(., outcome = "natdeaths")

hivdeaths <- lapply(bwaggr, attr, "hivdeaths") %>%
  abind::abind(., rev.along=0, use.dnns=TRUE,
               new.names = c(dimnm, list(sampleid = seq_along(.)))) %>%
  melt %>%
  data.table(., outcome = "hivdeaths")

outputs_age <- rbind(totpop, hivpop, infections, natdeaths, hivdeaths)

outputs_age$outcome <- factor(outputs_age$outcome, c("totpop", "hivpop", "infections", "natdeaths", "hivdeaths"))

dcast(outputs_age[sampleid == 1], year + sex + age ~ outcome) %>%
  write.csv("bwaggr_ageoutputs_single-resample.csv", row.names=FALSE)

dcast(outputs_age[ , .(value = median(value)), .(year, sex, age, outcome)],
      year + sex + age ~ outcome) %>%
  write.csv("bwaggr_ageoutputs_median.csv", row.names=FALSE)

dcast(outputs_age[ , .(value = mean(value)), .(year, sex, age, outcome)],
      year + sex + age ~ outcome) %>%
  write.csv("bwaggr_ageoutputs_mean.csv", row.names=FALSE)
```

For comparisons, the following annual aggregate outputs are generated:

* HIV prevalence among age 15-49
* HIV incidence rate among age 15-49
* ART coverage among age 15+


```r

hivprev <- sapply(bwaggr, prev) %>%
  "dimnames<-"(list(year = 1970:2021, sampleid = seq_len(ncol(.)))) %>%
  melt %>%
  data.table(., outcome = "prev15to49")

hivincid <- sapply(bwaggr, incid) %>%
  "dimnames<-"(list(year = 1970:2021, sampleid = seq_len(ncol(.)))) %>%
  melt %>%
  data.table(., outcome = "incid15to49")

artcov <- sapply(bwaggr, artcov15plus) %>%
  "dimnames<-"(list(year = 1970:2021, sampleid = seq_len(ncol(.)))) %>%
  melt %>%
  data.table(., outcome = "artcov15plus")

outputs_aggr <- rbind(hivprev, hivincid, artcov)

outputs_aggr$outcome <- factor(outputs_aggr$outcome,
                               c("prev15to49", "incid15to49", "artcov15plus"))

dcast(outputs_aggr[sampleid == 1], year ~ outcome) %>%
  write.csv("bwaggr_aggroutputs_single-resample.csv", row.names=FALSE)

dcast(outputs_aggr[ , .(value = median(value)), .(year, outcome)],
      year ~ outcome) %>%
  write.csv("bwaggr_aggroutputs_median.csv", row.names=FALSE)

dcast(outputs_aggr[ , .(value = mean(value)), .(year, outcome)],
      year ~ outcome) %>%
  write.csv("bwaggr_aggroutputs_mean.csv", row.names=FALSE)
```

## Parameter vectors for single resample results

Storing the paremter vectors for the single resample outputs so that we can 
reproduce and step through the simulations if necessary.


```r
theta_urban1 <- dput(c(bwfit$Urban$resample[1,]))
#> c(-0.811124446552004, -2.89213496664821, -1.25235012513843, 1997.27497455633, 
#> 0.135162964320709, 0.0897769226332093, 0.0766672488074031, -0.0390064645413289, 
#> -0.0276654028058017, -0.0684994868067683, 0.0706483722500527, 
#> 0.0123485242938665, -0.0500008566848367, 0.0242157988631865, 
#> -0.0213526437264166, 0.13014491136253, -0.0878701619022361, -0.0106679054280976, 
#> 0.0136766006923482, -0.00528706253043161, 0.0289274297466378, 
#> 1.55644249792549, -0.0753853509955008, -4.81924997735488, 0.00953586913866561, 
#> -7.40506288437238, -0.0191409304750598, 0.632085821778265, -1.84722193368205, 
#> 0.301224951856191, 0.613658351208605, 0.512955448252942, -0.0517366257251735, 
#> -0.531228949925646, 0.28769604878609, 0.213516323738396, -0.553827666537116, 
#> -0.558655731386, -0.488487410503282, -0.809969311740722, -0.0732092817360879, 
#> -0.0233701454601598, -0.151473583489263, -0.021302681885315, 
#> -0.00273859494434153, -0.0641227156782755)

theta_rural1 <- dput(c(bwfit$Rural$resample[1,]))
#> c(-0.970117400155276, -2.93469673090317, -0.66889538823997, 1999.62966842386, 
#> 0.0647120344521509, 0.00791627327988936, -0.0459624093224883, 
#> 0.014908161411549, -0.00668542243650456, 0.00411390841945129, 
#> 0.0203768672996544, 0.0298099710573834, 0.0209046894087502, 0.0048188871940908, 
#> 0.0504480205246233, -0.0328678900679028, -0.128017930978359, 
#> -0.0780009281325646, -0.031744274678965, 0.0460167468814374, 
#> -0.0712721596895912, 1.65068191620971, 0.0189361760294208, -5.24887127362293, 
#> 0.462604244503283, -9.90840310276572, 0.0114397423024262, 0.398753052034195, 
#> -1.76809122150709, -0.230444827099253, 0.347952253331402, 0.292622060217203, 
#> -0.132743218529816, -0.361746789893855, -0.220442316101592, 0.15342462060276, 
#> -0.176698024607871, -0.343068342736952, -0.990890966615118, -0.83661379715799, 
#> 0.121247572788112, 0.0696620905918071, -0.0120576304075231, -0.0940229336872126, 
#> 0.0477155054883171, -0.0572580437117387)
```
