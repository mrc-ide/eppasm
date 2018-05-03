#' ---
#' title: "Use of EPP-ASM for CSAVR"
#' author: Jeff Eaton, Guy Mahiane, Rob Glaubius, John Stover
#' output: pdf_document
#' ---
#'

##+ setup, include=FALSE
library(knitr)
opts_chunk$set(tidy=TRUE, warning=FALSE, cache=TRUE, message=FALSE)
options(knitr.kable.NA = '')

##+ load packages, include=FALSE
devtools::install_github("mrc-ide/epp",auth_token = pat_token)
devtools::install_github("mrc-ide/eppasm@csavr",auth_token = pat_token)
githubinstall::gh_install_packages("mrc-ide/eppasm",ref = "csavr",auth_token=pat_token)
library(eppasm)
## devtools::load_all("~/Documents/Code/R/eppasm-csavr/") # @csavr
devtools::load_all("hiv_project/jeff_eppasm_csavr_branch/")
devtools::build("hiv_project/jeff_eppasm_csavr_branch/")

library(epp)


library(magrittr)
library(broom)
library(ggplot2)


#' # Introduction
#'
#' Following discussion at the October 2017 UNAIDS Reference Group meeting, we showed
#' that the code for the EPP-ASM model was able to take annual HIV incidence rates as
#' inputs and reproduce Spectrum model ouputs for the numbers of new HIV infections,
#' PLHIV, and AIDS deaths for the age 15+ population. This report extends that
#' analysis to explore combining the EPP-ASM and CSAVR models for epidemic inference
#' from case surveillance and vital registration data.
#'
#' There a number of motivations for exploring a combined model:
#' 
#' * Consider opportunities for more efficient estimation for CSAVR countries.
#' * Explore influence of modelling incidence rate or transmission rate for
#'   estimates of recent incidence trends and potential harmonization of assumptions
#'   for CSAVR and EPP.
#' * Work towards incorporating case surveillance data into HIV estimates for SSA.
#' 
#' # Example data: Netherlands and Chile 2017 Spectrum files
#'
#' As an example case study, we use data about new diagnoses and numbers of
#' AIDS deaths from the Netherlands and Chile 2017 Spectrum files. For the Netherlands,
#' the Spectrum file only includes a partial time series of AIDS deaths, so we updated
#' these with the full time series of AIDS deaths to adults aged 15+ from the WHO
#' Mortality database. Moreover, Netherlands uses the ECDC model for direct incidence
#' estimates, not CSAVR, and so it is uncertain the source or accuracy of new case
#' inputs in the Spectrum file. In short, consider any outputs a proof of concept only.
#'

##+ read data, include=FALSE
## Read Spectrum inputs for direct incidence simulation
nl_pjnz <- "hiv_project/jeff_eppasm_data/Netherlands_2017_final.PJNZ"

nl_fp <- prepare_directincid(nl_pjnz)
nl_mod <- simmod(nl_fp)

nl_fp$relinfectART <- 0.3
nl_fp$tsEpidemicStart <- 1970.5

cl_pjnz <- "hiv_project/jeff_eppasm_data/Chile_2017_final.pjnz"

cl_fp <- prepare_directincid(cl_pjnz)
cl_mod <- simmod(cl_fp)

cl_fp$relinfectART <- 0.3
cl_fp$tsEpidemicStart <- 1970.5

## Read CSAVR data 
nl_csavrd <- read_csavr_data(nl_pjnz)
nl_csavrd[nl_csavrd == 0] <- NA
nl_csavrd$idx <- nl_csavrd$year - nl_fp$ss$proj_start + 1L

## Update Netherlands data with AIDS deaths from WHO mortality database
nl_who_aidsdeaths <- setNames(c(8, 16, 30, 63, 105, 132, 201, 267, 291, 409, 424, 444, 435, 320, 181, 134,
                                137, 132, 127, 89, 87, 85, 80, 48, 66, 53, 70, 50, 54, 44, 34, 38, 33),
                              1983:2015)
nl_csavrd[names(nl_who_aidsdeaths),]$aids_deaths <- nl_who_aidsdeaths

cl_csavrd <- read_csavr_data(cl_pjnz)
cl_csavrd[cl_csavrd == 0] <- NA
cl_csavrd$idx <- cl_csavrd$year - cl_fp$ss$proj_start + 1L

#' The tables below summarize the data for the Netherlands and Chile datasets.
##+ data summary, echo=FALSE
nl_csavrd[,c("plhiv", "plhiv_undercount", "new_cases", "new_cases_undercount", "aids_deaths", "aids_deaths_undercount")] %>%
  subset(!is.na(plhiv) | !is.na(new_cases) | !is.na(aids_deaths)) %>%
  setNames(sub(".*\\_undercount", "undercount", names(.))) %>%
  kable(caption = "Netherlands case surveillance and vital registration data")

cl_csavrd[,c("plhiv", "plhiv_undercount", "new_cases", "new_cases_undercount", "aids_deaths", "aids_deaths_undercount")] %>%
  subset(!is.na(plhiv) | !is.na(new_cases) | !is.na(aids_deaths)) %>%
  setNames(sub(".*\\_undercount", "undercount", names(.))) %>%
  kable(caption = "Chile case surveillance and vital registration data")


#' # Model for case surveillance data
#'
#' ## Model for new diagnoses
#'
#' For demonstration purposes, we used an much simplified model for new diagnoses.
#' We assumed that the untreated population not on ART is equivalent to the
#' undiagnosed population and modelled the diagnosis rate $\Delta_{s,a,m}(t)$ for
#' untreated HIV positive persons of sex $s$, age group $a$, and in CD4 stage $m$
#' at time $t$. The relative diagnosis rate across groups at time $t$ is assumed
#' to be proportional to the HIV mortality rate $\mu_{s,a,m}$ and the trend in
#' diagnosis rate of the course of the epidemic is modelled by a cumulative gamma
#' distribution function with shape parameter 1 and rate parameter $\theta$. That is
#' \begin{equation}
#'   \Delta_{s,a,m}(t) = \mu_{s,a,m} \cdot \gamma_{max} \cdot
#'   \int_{t_0}^t e^{-\theta \tau} d\tau
#' \end{equation}
#' where $t_0$ is the start time
#' of the epidemic (e.g. $t_0 = 1970$).
#'
#' Unlike the approximation by Mahiane, this simple model equating the untreated
#' population to the undiagnosed population does not track the proportion diagnosed
#' over the course of the epidemic. The full model should be implemented if this
#' approach is taken beyond this proof of concept.
#'
#' For Bayesian inference, we define diffuse prior distributions on the
#' parameters for the diagnosis rate over time:
#' $$ \log(\gamma_{max}) \sim \mathrm{normal}(3, 5) $$
#' $$ \log(\theta) \sim \mathrm{normal}(-3, 5) $$
#' 
#'
#' ## Models for incidence and transmission rate
#'
#' We considered three parameteric models for the HIV incidence rate, denoted
#' $\lambda(t)$, or HIV transmission $r(t)$.  For directly modelling $\lambda(t)$,
#' we implemented the single logistic (two parameters) and double logistic (five
#' parameters) models implemented by the current CSAVR software. The single logistic
#' model is
#' \begin{equation}
#'   \lambda(t) = c \cdot \frac{e^{\alpha(t-t_0)}}{1+e^{\alpha(t-t_0)}}
#' \end{equation}
#' where $c$ is the equilibrium incidence rate and $alpha$ determines the rate of
#' increase in incidence. The double logistic model is
#' \begin{equation}
#'    \lambda(t) = \frac{e^{\alpha(t-t_{mid})}}{1 + {e^{\alpha(t-t_{mid})}}} \cdot
#'    \left(2\cdot a \cdot \frac{e^{\beta(t-t_{mid})}}{1 + {e^{\beta(t-t_{mid})}}}
#'    + b\right)
#' \end{equation}
#' $b$ is the equilibrium incidence reate, $(a+b)/2$ is the incidence rate at the
#' inflection point $t_{mid}$, $\alpha$ determines the rate of increase and $\beta$
#' determines the rate of convergence to the asymptote.
#'
#' Secondly, instead of directly modelling the HIV incidence rate $\lambda(t)$,
#' we considered modelling the transmission rate $r(t)$, as in the EPP model. In
#' this case, the incidence rate is 
#' \begin{equation}
#'   \lambda(t) = r(t) \cdot \frac{I(t)}{N(t)} \cdot \left(1 - 0.7 * \frac{A(t)}{I(t)}\right).
#' \end{equation}
#' The expression $I(t)/N(t)$ is the HIV prevalence at time $t$, $A(t)/I(t)$ is the
#' ART coverage and 0.7 is the average reduction in transmission per additional person
#' on ART. We use a logistic function to model the logarithm or $r(t)$, termed the
#' *rlogistic* model, with four parameters
#' \begin{equation}
#'   \log r(t) = r_0 - (r_{\infty} - r_0)\cdot\frac{1}{1 + e^{-\alpha \cdot (t - t_{mid})}}
#' \end{equation}
#' where $e^{r_0}$ is the initial exponential growth rate of the epidemic,
#' $e^{r_{\infty}}$ is the equilibrium value for $r(t)$, $\alpha$ is the rate
#' of change in $\log r(t)$ and $t_{mid}$ is the inflection point. For this model
#' we additionally specify a fifth parmaeter $\iota$ as the incidence rate at time
#' $t=t_0$ providing the initial pulse of infections.
#' 
#' For Bayesian inference, we defined diffuse prior distributions on all parameters.
#' For the single and double logistic models, we defined the following prior 
#' distributions
#' $$ \log \alpha, \log \beta \sim \mathrm{normal}(-1, 5) $$
#' $$ \log a, \log b \sim \mathrm{normal}(-10, 5) $$
#' $$ t_{mid} \sim \mathrm{normal}(1995, 10) $$
#'
#' For the rlogistic model, we used prior distributions
#' $$ r_0 \sim \mathrm{normal}(\log(0.35), 0.5) $$
#' $$ r_\infty \sim \mathrm{normal}(\log(0.09), 0.3) $$
#' $$ \log \alpha \sim \mathrm{normal}(\log(0.2), 0.5) $$
#' $$ t_{mid} \sim \mathrm{normal}(1993, 5) $$
#' $$ \log\iota \sim \mathrm{normal}(-13, 5) $$ 
#'
#' Semi-parametric model variants for incidence rate or transmission rate (segmented
#' polynomials, p-splines, random walk, etc.) remain to be implemented and tested.
#' 
#' ## Likelihood
#'
#' The total number of expected diagnoses in year $t$ is the sum of expected diagnoses
#' over all sex, age, and CD4 groups of undiagnosed persons
#' $$ n_I(t) = \sum_s\sum_a\sum_m \Delta_{s,a,m}(t) \cdot U_{s,a,m}(t) $$
#' where $\Delta_{s,a,m}(t)$ is the diagnosis rate describe above and $U_{s,a,m}(t)$
#' are the number HIV positive undiagnosed, recalling that for simplicity we assumed
#' this was equal to the untreated population. A Poisson distribution is used as the
#' likelihood for the reported number of diagnoses $y_I(t)$ given the expected number
#' of diagnoses $n_I(t)$ and the proportion estimated underreporting of new diagnoses
#' $u_I(t)$
#' \begin{equation}
#'   y_I(t) \sim \mathrm{Poisson}(n_I(t) \cdot (1 - u_I(t)))
#' \end{equation}
#'
#' A Poisson likelihood is also used for the reported number of AIDS deaths $y_D(t)$
#' given the expected number of AIDS deaths $n_D(t)$ predicted by the Spectrum model
#' and the estimated underreporting of AIDS deaths $u_D(t)$:
#' \begin{equation}
#'   y_D(t) \sim \mathrm{Poisson}(n_D(t) \cdot (1 - u_D(t)))
#' \end{equation}
#'
#' We have not yet implemented the likelihood for the mean CD4 count at diagnosis
#' developed by Mahiane.
#'
#' # Model estimation
#'
#' To fit the model, we first create a model fitting object consisting of the
#' Spectrum model inputs and the CSAVR data inputs.
##+ prepare fitting
## Create fitting objects
nl <- list(fp = nl_fp, csavrd = nl_csavrd)
cl <- list(fp = cl_fp, csavrd = cl_csavrd)


#'
#' The code below illustrates fitting each of the three models to the Netherlands
#' dataset using either optimization (`optfit = TRUE`) or full Bayesian inference
#' with IMIS. The function `fitmod_csavr(...)` will prepare the model fit and data,
#' and then call the requested model fitting routine.

##+ fit nl, results = FALSE
## fit single logistic model for incidence rate
nl_opt1 <- fitmod_csavr(nl, incid_func = "ilogistic", B0=1e4, optfit=TRUE)
nl_fit1 <- fitmod_csavr(nl, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit double logistic model for incidence rate
nl_opt2 <- fitmod_csavr(nl, incid_func = "idbllogistic", B0=1e3, optfit=TRUE)
nl_fit2 <- fitmod_csavr(nl, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
nl_opt3 <- fitmod_csavr(nl, eppmod="rlogistic", B0=1e4, optfit=TRUE)
nl_fit3 <- fitmod_csavr(nl, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

#' Fit model to Chile dataset.

##+ fit cl, results = FALSE
## fit single logistic model for incidence rate
cl_opt1 <- fitmod_csavr(cl, incid_func = "ilogistic", B0=1e4, optfit=TRUE)
cl_fit1 <- fitmod_csavr(cl, incid_func = "ilogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit double logistic model for incidence rate
cl_opt2 <- fitmod_csavr(cl, incid_func = "idbllogistic", B0=1e3, optfit=TRUE)
cl_fit2 <- fitmod_csavr(cl, incid_func = "idbllogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)

## fit logistic model for transimssion rate (r(t))
cl_opt3 <- fitmod_csavr(cl, eppmod="rlogistic", B0=1e4, optfit=TRUE)
cl_fit3 <- fitmod_csavr(cl, eppmod="rlogistic", B0=1e4, B=1e3, B.re=3e3, opt_iter=1:3*5)


##+ tidy results, include = FALSE
nl_out1 <- tidy(nl_fit1) %>% data.frame(model = "logistic", .)
nl_out2 <- tidy(nl_fit2) %>% data.frame(model = "double logistic", .)
nl_out3 <- tidy(nl_fit3) %>% data.frame(model = "rlogistic", .)
nl_out <- rbind(nl_out1, nl_out2, nl_out3)

cl_out1 <- tidy(cl_fit1) %>% data.frame(model = "logistic", .)
cl_out2 <- tidy(cl_fit2) %>% data.frame(model = "double logistic", .)
cl_out3 <- tidy(cl_fit3) %>% data.frame(model = "rlogistic", .)
cl_out <- rbind(cl_out1, cl_out2, cl_out3)

#'
#' The figure below shows model results for the Netherlands from IMIS.
##+ plot nl fit, fig.height = 4.5, fig.width = 7, fig.align = "center", echo=FALSE
ggplot(subset(nl_out, year %in% 1975:2017), aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line() + geom_ribbon(linetype = "blank", alpha=0.25) + 
  facet_wrap(~outcome, scales="free") + 
  geom_point(aes(y=lik_data), col="darkred", size=1) +
  geom_point(aes(y=vld_data), col="grey40", size=1)+ 
  theme(legend.position = "bottom") +
  scale_x_continuous(element_blank()) + scale_y_continuous(element_blank())

#'
#' The figure below shows model results for Chile.
##+ plot cl fit, fig.height = 4.5, fig.width = 7, fig.align = "center", echo=FALSE 
ggplot(subset(cl_out, year %in% 1975:2017), aes(year, mean, ymin=lower, ymax=upper, color=model, fill=model)) +
  geom_line() + geom_ribbon(linetype = "blank", alpha=0.2) + 
  facet_wrap(~outcome, scales="free") + 
  geom_point(aes(y=lik_data), col="darkred", size=1) +
  geom_point(aes(y=vld_data), col="grey40", size=1)+ 
  theme(legend.position = "bottom") +
  scale_x_continuous(element_blank()) + scale_y_continuous(element_blank())

#' We should be cautious about over interpreting these results because we have not
#' fully implemented the model for new diagnoses and all of the parametric
#' models implemented thus far are relatively inflexible to matching data-driven
#' trends. A few observations:
#' 
#' * Both the double logistic model and rlogistic model, which each have five parameters,
#'   are better able to fit the data compared to the two parameter single logistic model.
#' * The uncertainty ranges are very narrow. This is likely due to both using relatively
#'   few parameters and not considering any uncertainty in other model processes such
#'   as disease progression or changes in testing and new diagnoses.
#' * Both the double logistic model and rlogistic model give fairly similar fits to
#'   observed new diagnoses and AIDS deaths, but result in qualitatively different
#'   estimates for recent HIV incidence trends. This conclusion should be reconsidered
#'   after implementing more flexible models for the incidence rate and transmission rate.
#'   
#' # Benchmarking
#'
#' The table summarizes the time in seconds for model fitting via optimization
#' and IMIS.
#'
#' 

##+ timings, echo=FALSE

kable(data.frame(Country = rep(c("Netherlands", "Chile"), each=3),
                 Model = rep(c("logistic", "double logistic", "r(t) logistic"), 2),
                 Optimization = c(nl_opt1$time[1], nl_opt2$time[1], nl_opt3$time[1],
                                  cl_opt1$time[1], cl_opt2$time[1], cl_opt3$time[1]),
                 IMIS = c(nl_fit1$time[1], nl_fit2$time[1], nl_fit3$time[1],
                          cl_fit1$time[1], cl_fit2$time[1], cl_fit3$time[1])),
      digits=c(NA, NA, 2, 1), align=c("llcc"))
                            
