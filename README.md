
<!-- README.md is generated from README.Rmd. Please edit that file -->

# serp <img src='man/figures/hex_logo.png' align="right" height="105" />

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being
activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://codecov.io/gh/ejikeugba/serp?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/serp)](https://CRAN.R-project.org/package=serp)
[![R build
status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)

<!-- badges: end -->

Smooth Effects on Response Penalty for CLM

A regularization method for the cumulative link models (CLM). The ‘serp’
function applies the ‘smooth-effect-on-response penalty’ (SERP) on the
estimates of the general CLM, enabling all subject-specific effects
associated with each variable in the model to shrink towards a unique
global effect. Fitting is done using a modified Newton’s method (Ugba et
al., 2021; Tutz and Gertheiss, 2016).

## Example

Consider the cumulative logit model of the wine dataset, where the
rating of wine bitterness is predicted with the two treatment factors,
temperature and contact.

``` r
require(serp)
#> Loading required package: serp
```

  - The unpenalized non-proportional odds model returns unbounded
    parameter estimates.

<!-- end list -->

``` r
require(serp)
m1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = TRUE, link = "logit", data = wine)
summary(m1)
#> 
#> call:
#> serp(formula = rating ~ temp + contact, link = "logit", slope = "unparallel", 
#>     reverse = TRUE, data = wine)
#> 
#> Coefficients:
#>                Estimate Std Error z value Pr(>|z|)    
#> (Intercept):1    1.2258    0.5570   2.201 0.027757 *  
#> (Intercept):2   -1.0330    0.4807  -2.149 0.031638 *  
#> (Intercept):3   -3.9464    0.9015  -4.377  1.2e-05 ***
#> (Intercept):4  -19.1843 1079.7764  -0.018 0.985825    
#> tempwarm:1      19.2427 3597.5553   0.005 0.995732    
#> tempwarm:2       2.1108    0.6007   3.514 0.000441 ***
#> tempwarm:3       2.9404    0.8283   3.550 0.000385 ***
#> tempwarm:4      17.0636 1079.7763   0.016 0.987392    
#> contactyes:1     1.6594    1.1786   1.408 0.159167    
#> contactyes:2     1.3429    0.5830   2.303 0.021258 *  
#> contactyes:3     1.6928    0.6602   2.564 0.010347 *  
#> contactyes:4     1.1622    0.9054   1.284 0.199293    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Number of iterations: 18 
#> 
#> Loglik: -84.61093 on 276 degrees of freedom 
#> 
#> AIC: 193.2219
#> 
#> Exponentiated coefficients:
#>   tempwarm:1   tempwarm:2   tempwarm:3   tempwarm:4 contactyes:1 contactyes:2 
#> 2.275154e+08 8.255243e+00 1.892362e+01 2.574225e+07 5.256007e+00 3.830073e+00 
#> contactyes:3 contactyes:4 
#> 5.434719e+00 3.196963e+00
```

  - Using SERP with the deviance tuning, for instance, returns the model
    along parameter shrinkage at which the total residual deviance is
    minimal and stable parameter estimates too.

<!-- end list -->

``` r

m2 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "deviance",
           data = wine)
summary(m2)
#> 
#> call:
#> serp(formula = rating ~ temp + contact, link = "logit", slope = "penalize", 
#>     tuneMethod = "deviance", reverse = TRUE, data = wine)
#> 
#> Coefficients:
#>               Estimate Std Error z value Pr(>|z|)    
#> (Intercept):1   1.2551    0.5558   2.258 0.023929 *  
#> (Intercept):2  -1.0613    0.4798  -2.212 0.026958 *  
#> (Intercept):3  -3.9263    0.8813  -4.455 8.39e-06 ***
#> (Intercept):4  -6.4972    2.4862  -2.613 0.008968 ** 
#> tempwarm:1      3.9218    2.3329   1.681 0.092751 .  
#> tempwarm:2      2.1494    0.5997   3.584 0.000338 ***
#> tempwarm:3      2.9272    0.8099   3.614 0.000301 ***
#> tempwarm:4      4.3677    2.4152   1.808 0.070538 .  
#> contactyes:1    1.6237    1.1257   1.442 0.149189    
#> contactyes:2    1.3675    0.5801   2.357 0.018404 *  
#> contactyes:3    1.6894    0.6519   2.592 0.009553 ** 
#> contactyes:4    1.1580    0.8922   1.298 0.194288    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Number of iterations: 10 
#> 
#> Loglik: -84.39441 on 276 degrees of freedom 
#> 
#> AIC: 192.7888
#> 
#> Exponentiated coefficients:
#>   tempwarm:1   tempwarm:2   tempwarm:3   tempwarm:4 contactyes:1 contactyes:2 
#>    50.489269     8.580056    18.674798    78.859536     5.071594     3.925343 
#> contactyes:3 contactyes:4 
#>     5.416215     3.183674 
#> 
#> Regularization Info:
#> penalty:     SERP
#> tuneMethod:  deviance
#> value:       168.7888
#> lambda:      0.07197
```

  - A user-supplied discrete lambda grid could be used for the deviance
    or the cv tuning method.

<!-- end list -->

``` r
m3 = serp(rating ~ temp + contact, slope = "penalize",
          reverse = TRUE, link = "logit", tuneMethod = "deviance",
          lambdaGrid = 10^seq(10, -2, length.out=10), data = wine)
head(predict(m3, type='response'))
#>            1         2         3          4           5
#> 1 0.21732519 0.5315983 0.2307966 0.01758655 0.002693383
#> 2 0.21732519 0.5315983 0.2307966 0.01758655 0.002693383
#> 3 0.05373011 0.3710352 0.4759543 0.09064459 0.008635883
#> 4 0.05373011 0.3710352 0.4759543 0.09064459 0.008635883
#> 5 0.01014334 0.2392808 0.4796472 0.16657436 0.104354292
#> 6 0.01014334 0.2392808 0.4796472 0.16657436 0.104354292
```

  - A penalized partial proportional odds model is obtained by setting
    one or variables as global effect(s).

<!-- end list -->

``` r
m4 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = TRUE, link = "logit", tuneMethod = "deviance",
           globalEff = ~ temp, data = wine)
coef(m4)
#> (Intercept):1 (Intercept):2 (Intercept):3 (Intercept):4      tempwarm 
#>      1.336395     -1.218164     -3.527260     -4.839052      2.414606 
#>  contactyes:1  contactyes:2  contactyes:3  contactyes:4 
#>      1.628993      1.526065      1.729784      1.418281
```

  - The Figure below shows the coefficient paths (row1) and the
    smoothing steps (row2) for the ordinal model of the wine data using
    SERP. The dashed horizontal (blue) and vertical (red) lines denote
    the parallel estimates and the selected estimates based on the
    minimum deviance, respectively.

<img src='docs/reference/figures/serpshrink.png' width="500" />

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

## References:

Ugba, E R; Mörlein, D; Gertheiss, J (2021). Smoothing in Ordinal
Regression: An Application to Sensory Data. *Stats*, 4, 616–633.
<https://doi.org/10.3390/stats4030037>

Tutz, G and Gertheiss, J (2016). Regularized Regression for Categorical
Data (With Discussion and Rejoinder). *Statistical Modelling*, 16,
161-260. <https://doi.org/10.1177/1471082X16642560>
