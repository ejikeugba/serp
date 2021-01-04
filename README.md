# serp
Smooth Effects on Response Penalty for CLM

A regularization method for the cumulative link models (CLM). 
The serp function applies the 'smooth-effect-on-response penalty' 
(SERP) on the estimates of the general CLM, enabling all 
subject-specific effects associated with each variable in 
the model to shrink towards a unique global effect.

## Example
Consider the cumulative logit model of the wine dataset,
where the rating of wine bitterness is predicted with 
the two treatment factors, temperature and contact. 

```R
## The unpenalized non-proportional odds model returns unbounded
## parameter estimates.
m1 <- serp(rating ~ temp + contact, slope = "unparallel",
           reverse = T, link = "logit", data = wine)
summary(m1)
 

## Using SERP with the deviance tuning,for instance,  returns 
## the model along parameter shrinkage at which the total 
## residual deviance is minimal and stable parameter estimates 
## too.
m2 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = T, link = "logit", tuning = "deviance",
           data = wine)
predict(m2)

## A user-supplied discrete lambda grid could be alternatively
## used for the deviance and cv methods.
m3 = serp(rating ~ temp + contact, slope = "penalize",
          reverse = T, link = "logit", tuning = "deviance",
          lambdagrid = seq(0,1,length.out = 10), data = wine)
confint(m3)


## A penalized partial proportional odds model is obtained by
## setting some variable(s) as global effect(s).
m4 <- serp(rating ~ temp + contact, slope = "penalize",
           reverse = T, link = "logit", tuning = "deviance",
           global = ~ temp, data = wine)
errorMetrics(m4)
```

## Installing:
```R
## install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

## Loading:
```R
library(serp)
