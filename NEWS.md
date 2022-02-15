## serp 0.2.3.9002
- errorMetrics() function is modified and no longer exported, an equivalent function with improved functionality is provided by the erroR() function in the gofcat package
- fix some dependency induced bugs in the package tests
- remove file License in the description and associated file docs
- minor changes in serp documentation

---
## serp 0.2.3
- CRAN release
- import from crayon, with colored outputs in returned objects 
- changes made in serp test for improved test coverage report
- bug fix in errorMetrics, with model argument also dropped 
- add print method and class of objects returned by errorMetrics
- message() replaces cat() where appropriate
- add citation for serp package

---
## serp 0.2.2
- JOSS release
- minor changes in serp documentation
- provide coefficient() as an alias for coef()
- update README.md to include community guidelines and contributors code of conduct
- 

---
## serp 0.2.1
- submission to CRAN

---
## serp 0.2.0.9001
- update function included in namespace
- examples included in the different function documentation
- shrinkage parameter upper limit set to 1e10 in serp.control
- bug fix in test function

---
## serp 0.2.0.9000 
- deviance tuning option in serp tuneMethod now replaced by AIC 
- serp output includes residual degrees of freedom (rdf)
- function to compute the effective degrees of freedom and rdf from the trace of the generalized hat matrix provided
- serp.summary documentation has its value segment edited
- changes made in serp test functions
- changes made in README.Rmd and README.md

---
## serp 0.2.0
- serp version 0.2.0 release

---
## serp 0.1.9.9001
- re-submission to CRAN

---
## serp 0.1.9.9000
- fixed all issues spotted out in the last version released on CRAN.  

---
## serp 0.1.9
- Submit serp version 0.1.9 to CRAN

---
## serp 0.1.8.9007
- updates README.md
- description gets additional doi
- references in serp documentation updated

---
## serp 0.1.8.9006
- Bugs in tests fixed

---
## serp 0.1.8.9005
- tests reconstructed to yield improved unit test coverage.
- rewrote some of the error and warning messages in key serp functions.
- Bugs in anova.serp fixed.
- 

---
## serp 0.1.8.9004
- update License in the description 
- Bugs in the example codes fixed
- Bugs in test files corrected

---
## serp 0.1.8.9003
- Column titles of predicted values with reverse TRUE are corrected.
- Bugs in the deviance (dvfun) are corrected.
- serp predict function now handles both single and multiple row input(s).
- serp.fit iteration algorithm improved.
- commented lines in serp main function removed.
- reverse and linkf arguments moved from main function to serpfit.
- loglog and cloglog links now give correct results with reverse='TRUE'.
- long lines of codes in the score function split for easy readability.
- starting values augmented to reflect the different link functions.
- reverse statement removed from cv function.
- reverse.fun is removed from summary.serp function.
- predicted values in errorMetrics get normalized for greater efficiency.

---
## serp 0.1.8.9002
- warning messages in serpfit get updated.
- Bug in serpfit when specifying gridType corrected.
- trainError no longer shows up in serp output.
- reverse.fun is now deprecated.
- penalty.print gets slightly reconstructed.
- Bugs in print.summary.serp fixed.
- Summary.serp drops off trainError with some few more adjustments.
- Predict gets reconstructed, also bugs when supplying 'newdata' corrected.
- Bugs in anova function fixed with slope type now included in the output.
- anova test for the 'penalize' slope currently disabled.
- updates in serp.control warning messages.

---
## serp 0.1.8.9001
- serp gets a new license.

---
## serp 0.1.8.9000
* README.md gets additional badges and a hexagon sticker

---
## serp 0.1.8
* package re-submission to CRAN

- reference to methods used in the package is included in the  description field
- values were added to all exported methods with corresponding explanations
- examples now runs by having all \dontrun{} removed
- all occurrences of <<- in the functions are dropped
- READmE.md is updated
