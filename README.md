
<!-- README.md is generated from README.Rmd. Please edit that file -->

# serp <img src='man/figures/hex_logo.png' align="right" height="105" />

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being
activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://codecov.io/gh/ejikeugba/serp?branch=master)
[![R build
status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/serp)](https://CRAN.R-project.org/package=serp)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](http://cranlogs.r-pkg.org/badges/grand-total)

<!-- badges: end -->

Smooth Effects on Response Penalty for CLM

A regularization method for the cumulative link models (CLM). The ‘serp’
function applies the ‘smooth-effect-on-response penalty’ (SERP) on the
estimates of the general CLM, enabling all subject-specific effects
associated with each variable in the model to shrink towards a unique
global effect.

## Example

Consider the cumulative logit model of the wine dataset, where the
rating of wine bitterness is predicted with the two treatment factors,
temperature and contact.

``` r
## The unpenalized non-proportional odds model returns unbounded
## parameter estimates.
m1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
summary(m1)
 

## Using SERP with the deviance tuning, for instance,  returns 
## the model along parameter shrinkage at which the total 
## residual deviance is minimal and stable parameter estimates 
## too.
m2 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "deviance",
           data = wine)
predict(m2, type='response')

## A user-supplied discrete lambda grid could be alternatively
## used for the deviance and cv methods.
m3 = serp(rating ~ temp + contact, slope = "penalize",
          reverse = TRUE, link = "logit", tuneMethod = "deviance",
          lambdaGrid = 10^seq(10, -2, length.out=10), data = wine)
confint(m3)


## A penalized partial proportional odds model is obtained by
## setting some variable(s) as global effect(s).
m4 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "deviance",
           globalEff = ~ temp, data = wine)
errorMetrics(m4)
```

## Installation:

You can install the released version of serp from
[CRAN](https://cran.r-project.org/package=serp) with:

``` r
install.packages("serp")
```

And the development version from
[GitHub](https://github.com/ejikeugba/serp) with:

``` r
# install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

## Loading:

``` r
library(serp)
```
