# Prepare design matrix indices for ANC prevalence predictions

Prepare design matrix indices for ANC prevalence predictions

## Usage

``` r
ancsite_pred_df(ancsite_df, fp)
```

## Arguments

- ancsite_df:

  data.frame of site-level ANC design for predictions

- fp:

  fixed parameter input list

## Examples

``` r
pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="eppasm")
bw <- prepare_spec_fit(pjnz, proj.end=2021.5)


bw_u_ancsite <- attr(bw$Urban, "eppd")$ancsitedat
fp <- attr(bw$Urban, "specfp")

ancsite_pred_df(bw_u_ancsite, fp)
#> $df
#>     agspan aidx yidx               site year used  prev    n  type agegr age
#> 1       35    1   22    Francistown (%) 1991 TRUE 0.160  128 ancss 15-49  15
#> 2       35    1   22       Gaborone (%) 1991 TRUE 0.170   58 ancss 15-49  15
#> 3       35    1   23       Gaborone (%) 1992 TRUE 0.149  841 ancss 15-49  15
#> 4       35    1   23    Francistown (%) 1992 TRUE 0.237  796 ancss 15-49  15
#> 5       35    1   24       Gaborone (%) 1993 TRUE 0.192  881 ancss 15-49  15
#> 6       35    1   24        Lobatse (%) 1993 TRUE 0.178  258 ancss 15-49  15
#> 7       35    1   24    Francistown (%) 1993 TRUE 0.342  803 ancss 15-49  15
#> 8       35    1   24 Serowe Palapye (%) 1993 TRUE 0.199  267 ancss 15-49  15
#> 9       35    1   25       Gaborone (%) 1994 TRUE 0.278 1205 ancss 15-49  15
#> 10      35    1   25       Southern (%) 1994 TRUE 0.160  508 ancss 15-49  15
#> 11      35    1   25    Francistown (%) 1994 TRUE 0.297  799 ancss 15-49  15
#> 12      35    1   25   Selebi-Pikwe (%) 1994 TRUE 0.270  307 ancss 15-49  15
#> 13      35    1   26        Lobatse (%) 1995 TRUE 0.389  231 ancss 15-49  15
#> 14      35    1   26    Francistown (%) 1995 TRUE 0.396  626 ancss 15-49  15
#> 15      35    1   26 Serowe Palapye (%) 1995 TRUE 0.299  262 ancss 15-49  15
#> 16      35    1   26       Gaborone (%) 1995 TRUE 0.287 1307 ancss 15-49  15
#> 17      35    1   27       Gaborone (%) 1996 TRUE 0.313 1232 ancss 15-49  15
#> 18      35    1   27   Kweneng East (%) 1996 TRUE 0.239  213 ancss 15-49  15
#> 19      35    1   27       Southern (%) 1996 TRUE 0.217  589 ancss 15-49  15
#> 20      35    1   27    Francistown (%) 1996 TRUE 0.431  751 ancss 15-49  15
#> 21      35    1   27   Selebi-Pikwe (%) 1996 TRUE 0.378  331 ancss 15-49  15
#> 22      35    1   28       Gaborone (%) 1997 TRUE 0.340  786 ancss 15-49  15
#> 23      35    1   28      Mahalapye (%) 1997 TRUE 0.282  393 ancss 15-49  15
#> 24      35    1   28        Lobatse (%) 1997 TRUE 0.337  266 ancss 15-49  15
#> 25      35    1   28    Francistown (%) 1997 TRUE 0.429  802 ancss 15-49  15
#> 26      35    1   28 Serowe Palapye (%) 1997 TRUE 0.344  289 ancss 15-49  15
#> 27      35    1   29   Kweneng East (%) 1998 TRUE 0.372  534 ancss 15-49  15
#> 28      35    1   29       Southern (%) 1998 TRUE 0.247  300 ancss 15-49  15
#> 29      35    1   29   Selebi-Pikwe (%) 1998 TRUE 0.499  471 ancss 15-49  15
#> 30      35    1   29    Francistown (%) 1998 TRUE 0.430  796 ancss 15-49  15
#> 31      35    1   29       Gaborone (%) 1998 TRUE 0.391 1251 ancss 15-49  15
#> 32      35    1   30       Gaborone (%) 1999 TRUE 0.371  480 ancss 15-49  15
#> 33      35    1   30        Lobatse (%) 1999 TRUE 0.313  252 ancss 15-49  15
#> 34      35    1   30    Francistown (%) 1999 TRUE 0.427  576 ancss 15-49  15
#> 35      35    1   30 Serowe Palapye (%) 1999 TRUE 0.418  268 ancss 15-49  15
#> 36      35    1   30      Mahalapye (%) 1999 TRUE 0.320  362 ancss 15-49  15
#> 37      35    1   31   Kweneng East (%) 2000 TRUE 0.304  349 ancss 15-49  15
#> 38      35    1   31   Selebi-Pikwe (%) 2000 TRUE 0.503  304 ancss 15-49  15
#> 39      35    1   31       Southern (%) 2000 TRUE 0.407  513 ancss 15-49  15
#> 40      35    1   31       Gaborone (%) 2000 TRUE 0.362  734 ancss 15-49  15
#> 41      35    1   31    Francistown (%) 2000 TRUE 0.444  702 ancss 15-49  15
#> 42      35    1   32 Serowe Palapye (%) 2001 TRUE 0.446  341 ancss 15-49  15
#> 43      35    1   32       Southern (%) 2001 TRUE 0.340  347 ancss 15-49  15
#> 44      35    1   32   Kweneng East (%) 2001 TRUE 0.296  540 ancss 15-49  15
#> 45      35    1   32     South East (%) 2001 TRUE 0.321  192 ancss 15-49  15
#> 46      35    1   32    Francistown (%) 2001 TRUE 0.496  494 ancss 15-49  15
#> 47      35    1   32        Lobatse (%) 2001 TRUE 0.306  232 ancss 15-49  15
#> 48      35    1   32   Selebi-Pikwe (%) 2001 TRUE 0.500  351 ancss 15-49  15
#> 49      35    1   32       Gaborone (%) 2001 TRUE 0.386  576 ancss 15-49  15
#> 50      35    1   32      Mahalapye (%) 2001 TRUE 0.319  336 ancss 15-49  15
#> 51      35    1   33       Southern (%) 2002 TRUE 0.331  363 ancss 15-49  15
#> 52      35    1   33   Kweneng East (%) 2002 TRUE 0.292  567 ancss 15-49  15
#> 53      35    1   33     South East (%) 2002 TRUE 0.265  288 ancss 15-49  15
#> 54      35    1   33   Selebi-Pikwe (%) 2002 TRUE 0.481  372 ancss 15-49  15
#> 55      35    1   33    Francistown (%) 2002 TRUE 0.402  458 ancss 15-49  15
#> 56      35    1   33       Gaborone (%) 2002 TRUE 0.382  666 ancss 15-49  15
#> 57      35    1   33        Lobatse (%) 2002 TRUE 0.346  241 ancss 15-49  15
#> 58      35    1   33 Serowe Palapye (%) 2002 TRUE 0.367  425 ancss 15-49  15
#> 59      35    1   33      Mahalapye (%) 2002 TRUE 0.398  303 ancss 15-49  15
#> 60      35    1   34   Selebi-Pikwe (%) 2003 TRUE 0.522  302 ancss 15-49  15
#> 61      35    1   34       Gaborone (%) 2003 TRUE 0.448  555 ancss 15-49  15
#> 62      35    1   34       Southern (%) 2003 TRUE 0.257  359 ancss 15-49  15
#> 63      35    1   34    Francistown (%) 2003 TRUE 0.458  593 ancss 15-49  15
#> 64      35    1   34 Serowe Palapye (%) 2003 TRUE 0.433  389 ancss 15-49  15
#> 65      35    1   34        Lobatse (%) 2003 TRUE 0.324  246 ancss 15-49  15
#> 66      35    1   34      Mahalapye (%) 2003 TRUE 0.374  326 ancss 15-49  15
#> 67      35    1   34     South East (%) 2003 TRUE 0.279  306 ancss 15-49  15
#> 68      35    1   34   Kweneng East (%) 2003 TRUE 0.321  522 ancss 15-49  15
#> 69      35    1   36      Mahalapye (%) 2005 TRUE 0.362  379 ancss 15-49  15
#> 70      35    1   36        Jwaneng (%) 2005 TRUE 0.349  146 ancss 15-49  15
#> 71      35    1   36   Kweneng East (%) 2005 TRUE 0.315  551 ancss 15-49  15
#> 72      35    1   36       Gaborone (%) 2005 TRUE 0.344  722 ancss 15-49  15
#> 73      35    1   36       Southern (%) 2005 TRUE 0.282  313 ancss 15-49  15
#> 74      35    1   36     South East (%) 2005 TRUE 0.321  368 ancss 15-49  15
#> 75      35    1   36        Lobatse (%) 2005 TRUE 0.290  243 ancss 15-49  15
#> 76      35    1   36    Francistown (%) 2005 TRUE 0.423  659 ancss 15-49  15
#> 77      35    1   36 Serowe Palapye (%) 2005 TRUE 0.375  534 ancss 15-49  15
#> 78      35    1   36   Selebi-Pikwe (%) 2005 TRUE 0.465  300 ancss 15-49  15
#> 79      35    1   37     South East (%) 2006 TRUE 0.276  385 ancss 15-49  15
#> 80      35    1   37   Kweneng East (%) 2006 TRUE 0.317  518 ancss 15-49  15
#> 81      35    1   37      Mahalapye (%) 2006 TRUE 0.299  322 ancss 15-49  15
#> 82      35    1   37        Jwaneng (%) 2006 TRUE 0.200  143 ancss 15-49  15
#> 83      35    1   37   Selebi-Pikwe (%) 2006 TRUE 0.411  324 ancss 15-49  15
#> 84      35    1   37       Southern (%) 2006 TRUE 0.250  324 ancss 15-49  15
#> 85      35    1   37        Lobatse (%) 2006 TRUE 0.284  249 ancss 15-49  15
#> 86      35    1   37       Gaborone (%) 2006 TRUE 0.353  607 ancss 15-49  15
#> 87      35    1   37 Serowe Palapye (%) 2006 TRUE 0.368  598 ancss 15-49  15
#> 88      35    1   37    Francistown (%) 2006 TRUE 0.367  560 ancss 15-49  15
#> 89      35    1   38   Selebi-Pikwe (%) 2007 TRUE 0.490  310 ancss 15-49  15
#> 90      35    1   38     South East (%) 2007 TRUE 0.267  441 ancss 15-49  15
#> 91      35    1   38        Jwaneng (%) 2007 TRUE 0.214  141 ancss 15-49  15
#> 92      35    1   38 Serowe Palapye (%) 2007 TRUE 0.373  667 ancss 15-49  15
#> 93      35    1   38       Southern (%) 2007 TRUE 0.228  319 ancss 15-49  15
#> 94      35    1   38   Kweneng East (%) 2007 TRUE 0.335  487 ancss 15-49  15
#> 95      35    1   38      Mahalapye (%) 2007 TRUE 0.304  403 ancss 15-49  15
#> 96      35    1   38       Gaborone (%) 2007 TRUE 0.332  692 ancss 15-49  15
#> 97      35    1   38        Lobatse (%) 2007 TRUE 0.287  242 ancss 15-49  15
#> 98      35    1   38    Francistown (%) 2007 TRUE 0.379  702 ancss 15-49  15
#> 99      35    1   40        Jwaneng (%) 2009 TRUE 0.311  300 ancss 15-49  15
#> 100     35    1   40   Kweneng East (%) 2009 TRUE 0.279  486 ancss 15-49  15
#> 101     35    1   40 Serowe Palapye (%) 2009 TRUE 0.336  405 ancss 15-49  15
#> 102     35    1   40   Selebi-Pikwe (%) 2009 TRUE 0.391  263 ancss 15-49  15
#> 103     35    1   40        Lobatse (%) 2009 TRUE 0.324  274 ancss 15-49  15
#> 104     35    1   40     South East (%) 2009 TRUE 0.288  421 ancss 15-49  15
#> 105     35    1   40    Francistown (%) 2009 TRUE 0.381  527 ancss 15-49  15
#> 106     35    1   40       Gaborone (%) 2009 TRUE 0.290  678 ancss 15-49  15
#> 107     35    1   40       Southern (%) 2009 TRUE 0.283  316 ancss 15-49  15
#> 108     35    1   40      Mahalapye (%) 2009 TRUE 0.378  409 ancss 15-49  15
#> 109     35    1   42     South East (%) 2011 TRUE 0.158  242 ancss 15-49  15
#> 110     35    1   42    Francistown (%) 2011 TRUE 0.370  547 ancss 15-49  15
#> 111     35    1   42        Jwaneng (%) 2011 TRUE 0.167  300 ancss 15-49  15
#> 112     35    1   42   Selebi-Pikwe (%) 2011 TRUE 0.412  282 ancss 15-49  15
#> 113     35    1   42   Kweneng East (%) 2011 TRUE 0.266  262 ancss 15-49  15
#> 114     35    1   42      Mahalapye (%) 2011 TRUE 0.335  463 ancss 15-49  15
#> 115     35    1   42 Serowe Palapye (%) 2011 TRUE 0.379  458 ancss 15-49  15
#> 116     35    1   42       Gaborone (%) 2011 TRUE 0.292  572 ancss 15-49  15
#> 117     35    1   42        Lobatse (%) 2011 TRUE 0.244  169 ancss 15-49  15
#> 118     35    1   42       Southern (%) 2011 TRUE 0.224  302 ancss 15-49  15
#>     qMidx
#> 1       1
#> 2       1
#> 3       2
#> 4       2
#> 5       3
#> 6       3
#> 7       3
#> 8       3
#> 9       4
#> 10      4
#> 11      4
#> 12      4
#> 13      5
#> 14      5
#> 15      5
#> 16      5
#> 17      6
#> 18      6
#> 19      6
#> 20      6
#> 21      6
#> 22      7
#> 23      7
#> 24      7
#> 25      7
#> 26      7
#> 27      8
#> 28      8
#> 29      8
#> 30      8
#> 31      8
#> 32      9
#> 33      9
#> 34      9
#> 35      9
#> 36      9
#> 37     10
#> 38     10
#> 39     10
#> 40     10
#> 41     10
#> 42     11
#> 43     11
#> 44     11
#> 45     11
#> 46     11
#> 47     11
#> 48     11
#> 49     11
#> 50     11
#> 51     12
#> 52     12
#> 53     12
#> 54     12
#> 55     12
#> 56     12
#> 57     12
#> 58     12
#> 59     12
#> 60     13
#> 61     13
#> 62     13
#> 63     13
#> 64     13
#> 65     13
#> 66     13
#> 67     13
#> 68     13
#> 69     14
#> 70     14
#> 71     14
#> 72     14
#> 73     14
#> 74     14
#> 75     14
#> 76     14
#> 77     14
#> 78     14
#> 79     15
#> 80     15
#> 81     15
#> 82     15
#> 83     15
#> 84     15
#> 85     15
#> 86     15
#> 87     15
#> 88     15
#> 89     16
#> 90     16
#> 91     16
#> 92     16
#> 93     16
#> 94     16
#> 95     16
#> 96     16
#> 97     16
#> 98     16
#> 99     17
#> 100    17
#> 101    17
#> 102    17
#> 103    17
#> 104    17
#> 105    17
#> 106    17
#> 107    17
#> 108    17
#> 109    18
#> 110    18
#> 111    18
#> 112    18
#> 113    18
#> 114    18
#> 115    18
#> 116    18
#> 117    18
#> 118    18
#> 
#> $datgrp
#>    aidx yidx agspan qMidx
#> 1     1   22     35     1
#> 2     1   23     35     2
#> 3     1   24     35     3
#> 4     1   25     35     4
#> 5     1   26     35     5
#> 6     1   27     35     6
#> 7     1   28     35     7
#> 8     1   29     35     8
#> 9     1   30     35     9
#> 10    1   31     35    10
#> 11    1   32     35    11
#> 12    1   33     35    12
#> 13    1   34     35    13
#> 14    1   36     35    14
#> 15    1   37     35    15
#> 16    1   38     35    16
#> 17    1   40     35    17
#> 18    1   42     35    18
#> 
```
