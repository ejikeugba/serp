
<!-- README.md is generated from README.Rmd. Please edit that file -->

# serp <a href="https://ejikeugba.github.io/serp/"><img src='man/figures/hex_logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being
activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ejikeugba/serp?branch=master)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](https://CRAN.R-project.org/package=serp)
[![CRAN
status](https://www.r-pkg.org/badges/version/serp)](https://CRAN.R-project.org/package=serp)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03705/status.svg)](https://doi.org/10.21105/joss.03705)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![R build
status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)

<!-- badges: end -->

### Overview

The `serp` R package fits cumulative link models (CLMs) with the
`smooth-effect-on-response penalty (SERP)`. The `cumulative model`
developed by McCullagh (1980) is probably the most frequently used
ordinal model in empirical studies. However, the stochastic ordering
property of the general form of the model poses a very serious challenge
in most empirical applications of the model. For instance, unstable
likelihoods with ill-conditioned parameter space are frequently
encountered during the iterative process. `serp` implements a unique
regularization method for CLMs that provides the means of smoothing the
adjacent categories in the model. At extreme shrinkage, SERP causes all
subject-specific effects associated with each variable in the model to
shrink towards unique global effects. Fitting is done using a modified
Newton’s method. Several standard model performance and descriptive
methods are also available. See [Ugba,
2021](https://www.researchgate.net/publication/355796737_serp_An_R_package_for_smoothing_in_ordinal_regression), [Ugba et al.,
2021](https://www.researchgate.net/publication/353406284_Smoothing_in_Ordinal_Regression_An_Application_to_Sensory_Data) and [Tutz and Gertheiss,
2016](https://doi.org/10.1177/1471082X16642560) for further details on
the implemented penalty.

### Example

Consider the cumulative logit model of the [wine
dataset](https://ejikeugba.github.io/serp/reference/wine.html), where
the rating of wine bitterness is predicted with the two treatment
factors, temperature and contact.

``` r
## The unpenalized non-proportional odds model returns unbounded estimates, hence,
## not fully identifiable.
f1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
coef(f1)
```

``` r
## The penalized non-proportional odds model with a user-supplied lambda gives 
## a fully identified model having bounded estimates. A suitable tuning criterion
## could as well be used to select lambda (e.g., aic or cv) 
f2 <- serp(rating ~ temp + contact, slope = "penalize",
           link = "logit", reverse = TRUE, tuneMethod = "user",
           lambda = 1e1 ,data = wine)
coef(f2)
```

``` r
## A penalized partial proportional odds model with one variable set to 
## global effect is also possible.
f3 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "user",
           lambda = 2e1, globalEff = ~ temp, data = wine)
coef(f3)
```

``` r
## The unpenalized proportional odds model with constrained estimates. 
## Under estreme shrinkage, estimates in f2 equal those in this model.  
f4 <-  serp(rating ~ temp + contact, slope = "parallel",
            reverse = FALSE, link = "logit", data = wine)
summary(f4)
```

### Installation and Use

Before installing `serp`, it is encouraged to have a recent version of
[R](https://cran.r-project.org/bin/windows/base/) installed. The
released version of `serp` can be installed from
[CRAN](https://cran.r-project.org/package=serp) with:

``` r
install.packages("serp")
```

or the development version from
[GitHub](https://github.com/ejikeugba/serp) with:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

Load `serp` into R environment with:

``` r
library(serp)
```

### Community Guidelines

Pull requests are welcomed! Please submit your contributions to `serp`
through the list of `Pull Requests`, following the [contributing
guidelines](https://ejikeugba.github.io/serp/CONTRIBUTING.html). To
report issues and/or seek support, please file a new ticket in the
[issue](https://github.com/ejikeugba/serp/issues) tracker, and expect a
feedback ASAP!

### Code of Conduct

Please note that `serp` is released with a [Contributor Code of
Conduct](https://github.com/ejikeugba/serp/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

### References

McCullagh, P. (1980). Regression Models for Ordinal Data. *Journal of
the Royal Statistical Society. Series B (Methodological)*, 42, 109-142.
<https://doi.org/10.1111/j.2517-6161.1980.tb01109.x>

Randall, J (1989). The analysis of sensory data by generalized linear
model. *Biometrical Journal*, 31, 781–793.
<https://doi.org/10.1002/bimj.4710310703>

Tutz, G. and Gertheiss, J. (2016). Regularized Regression for
Categorical Data (With Discussion and Rejoinder). *Statistical
Modelling*, 16, 161-260. <https://doi.org/10.1177/1471082X16642560>

Ugba, E. R., Mörlein, D. and Gertheiss, J. (2021). Smoothing in Ordinal
Regression: An Application to Sensory Data. *Stats*, 4, 616–633.
<https://doi.org/10.3390/stats4030037>

Ugba, E. R. (2021). serp: An R package for smoothing in ordinal
regression *Journal of Open Source Software*, 6(66), 3705.
<https://doi.org/10.21105/joss.03705>
