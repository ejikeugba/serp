# serp
Smooth Effects on Response Penalty for CLM

A regularization method for the cumulative link models (CLM). 
The serp function applies the 'smooth-effect-on-response penalty' 
(SERP) on the estimates of the general CLM, forcing all 
subject-specific effects associated with each response category
to shrink towards a unique global effect.

## Installing:
```R
## install.packages("devtools")
devtools::install_github("ejikeugba/serp")
```

## Loading:
```R
library(serp)
