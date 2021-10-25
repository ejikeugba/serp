---
title: 'serp: An R package for smoothing in ordinal regression'
tags:
  - R
  - regularization
  - identification problem
  - cumulative models
  - categorical data
  - proportional odds
  - shrinkage penalty
authors:
  - name: Ejike R. Ugba
    orcid: 0000-0003-2572-0023
    affiliation: 1
affiliations:
 - name: Department of Mathematics and Statistics, School of Economics and Social Sciences, Helmut Schmidt University, Hamburg, Germany
   index: 1
date: 25 August 2021
bibliography: paper.bib
---


# Summary 

The use of specialized methods for the analysis of categorical data has been on the rise in recent years [@agresti_analysis_2010]. For instance, scientists frequently use the different forms of the ordinal model to analyze relationships between an ordinal response variable and covariates of interest. In particular, the cumulative link model (CLM) propagated by @mcCullagh_regression_1980 finds a wide range of empirical applications in clinical trials, social surveys, market research, etc. However, in high-dimensional settings with a large number of unknown parameters, e.g., if several potential predictors and response categories are modeled, the so-called identification problem could be somewhat unavoidable [@fisher_identification_1966; @bartels_identification_1985]. Thankfully, regularization techniques [see, e.g., @bhlmann_statistics_2011;@hastie_2009] among other approaches, offer a remedy in such situations.    

The `serp` software package [@Ugba_serp_2021], available in the Comprehensive R Archive Network ([CRAN](https://CRAN.R-project.org/package=serp)) [@R_language_2021], implements a unique regularization algorithm that enforces smoothness on the adjacent categories of CLMs. The approach provides flexible modeling of CLMs that uses the smooth-effects-on-response penalty (SERP) discussed in @ugba_smoothing_2021 and @tutz_regularized_2016. For instance, when applied to the general form of the cumulative logit model, also known as the non-proportional odds model (NPOM), one obtains the proportional odds model (POM) under extreme parameter shrinkage. Fitting in `serp` is carried out by a modified Newton's optimization method that facilitates an easy convergence and estimation speed. As an open-source software package, `serp` aims to provide a platform for advanced scientific modeling in empirical research. The core features of `serp`, as well as details about usage, are provided in the software package [documentation](https://cran.r-project.org/web/packages/serp/serp.pdf).



# Statement of Need

Although a very useful tool in regression analysis, the development of regularization methods for categorical response models is still in its elementary stage. This explains why so many popular statistical software packages for ordinal regression do not currently implement such methods. The 'vglm' function in the `VGAM` R package [@yee_vector_1996; @yee_vgam_2010], for instance, fits the cumulative family of the ordered model, with several link functions, but with no provision for regularization. The same is also true for the 'clm' and 'polr' functions in the `ordinal` [@christensen_ordinal_2019] and `MASS` [@venables_Modern_2002] R packages, respectively. Other non-R functions for ordinal regression, such as the SAS `CATMOD` procedure [@SAS_Institute_2018] and the SPSS `PLUM` procedure [@IBM_spss_2021] also do not currently support regularization. `serp` fills this gap by providing a very unique and convenient means of regularizing estimates in ordinal models. The software package is intended for diverse forms of scientific applications. It is particularly suited for studies involving ordinal categorical outcomes. Users already conversant with existing R packages for ordinal regression will by no means find it difficult to use `serp`. The package, moreover, comes with several useful functions that support empirical research while being also equipped with standard model performance and descriptive methods.



# State of the Field
As previously noted, although a topic of intensive research in statistics for decades, regularization methods that are specifically designed for categorical data are relatively new [@ugba_smoothing_2021]. So far, exiting software packages implementing such techniques are mainly extensions of methods typically found in the high-dimensional linear or generalized linear model framework. Mostly, penalties in the tradition of lasso [@tibshirani_regression_1996], ridge [@hoerl_ridge_1970] and elastic net [@zou_regularization_2005] are implemented. The `glmnetcr` [@archer_glmnetcr_2014a] and `glmpathcr` [@archer_glmpathcr_2014b] R packages, for instance, fit continuation ratio models with the elastic net penalty, while the `rms` R package [@harrell_regressionl_2021] fits the cumulative logit model with quadratic (ridge regression) penalty. Also, the `ordinalgmifs` R package provided by @archer_ordinalgmifs_2014 implements a generalized monotone incremental forward stagewise (GMIFS) algorithm for regularized ordinal regression models, with a solution path akin to the L1 norm (lasso) penalty. Furthermore, the `ordinalNet` R package [@wurm_regularized_2021] fits the elementwise link multinomial-ordinal (ELMO) class of models with the elastic net using a coordinate descent algorithm.

A common principle behind the regularization procedures with the highlighted software packages is to make the model's maximum likelihood estimates biased a bit towards zero or fully shrunk to zero, of course, with the so-called bias-variance trade-off in mind [@Fortmann-Roe_2012]. However, the regularization method implemented in `serp` makes a paradigm shift from the traditional form of regularization where estimates are shrunk towards zero. `serp`  provides a unique form of penalization that shrinks the category effects in the non-proportional cumulative link model towards global effects. Under extreme shrinkage, this form of penalization results in the proportional odds model. A somewhat similar form of smoothing for categorical response data is provided in the `mgcv` R package [@wood_smoothing_2016] but for un-ordered categories and with smoothing done in the linear predictors. By and large, regularization with `serp` provides a means of reducing model complexity without necessarily chopping off important features in the model. 



# A Minimal Example
The [wine dataset](https://ejikeugba.github.io/serp/reference/wine.html) from @Randall_1989 available in `serp` software package represents the outcome of a factorial experiment on factors determining the bitterness of wine. Two treatment factors (temperature and contact) with two levels each are provided. Indeed, temperature (warm or cold) and contact (yes or no) between juice and skins can be controlled when crushing grapes during wine production. Nine judges each assessed wine from two bottles from each of the four treatment conditions, making a total of 72 observations. Rating of wine bitterness ranges from 1 = 'least bitter' to 5 = 'most bitter'. A typical modeling framework in this instance is an ordinal regression where the cumulative logit model is used to predict the scores of wine bitterness using the stated treatment factors (i.e., a simple case where all observations are aggregated over bottles and judges, compare e.g., the examples given in the `ordinal` R package manual [@Christensen_cumulative_2019]). However, the fitted NPOM (see, Table 1) using the `VGAM`-vglm function [@yee_vgam_2010], for instance,  is not fully identifiable. As observed, unbounded estimates with large standard errors ($\mbox{SE} > 10^3$) were obtained for covariate temperature with respect to the first and last response category, and also the last threshold (Intercept:4). In the actual sense, the absolute values of the estimates and the estimated standard errors for the unbounded parameters seem to diverge to $\infty$ as the stopping criteria of the iterative fitting procedure used become stricter and the number of iterations allowed increases [see, e.g., @Kosmidis_bias_2014]. Other software implementations of the same model without regularization run into a similar problem, forcing users to adopt a different modeling approach which may or may not be appropriate for the data. However, by incorporating SERP into the modeling framework, a fully identified model with bounded estimates could be achieved. The realized coefficient paths for the estimated model via `serp`, using increasing values of the tuning parameter $[\lambda: 10^{-3}, 10^5]$, are shown in \autoref{fig:serpfig}. 


![Estimated coefficient paths for the ordinal model of the wine data when using the smooth-effects-on-response penalty (SERP). The thick lines on the top displays are the category-specific coefficients associated with the two predictors, under increasing values of the tuning parameter $(\lambda)$ on the range ($10^{-3}, 10^{5}$). The dashed horizontal (blue) lines denote the parallel estimates. The bottom displays further illustrate SERP's smoothing steps from the category-specific to the parallel estimates, with the solid black, grey and dashed (blue) line strokes respectively the category-specific, the penalized and the parallel estimates. \label{fig:serpfig}](serp_fig.png){ width=70%}



  Coefficients  |        vglm        |        serp        |
  :-------------|-------------------:|-------------------:|
  (Intercept):1 |   1.226    (0.557) |    1.344   (0.509) |
  (Intercept):2 |  -1.033    (0.481) |   -1.251   (0.439) | 
  (Intercept):3 |  -3.946    (0.902) |   -3.467   (0.597) |
  (Intercept):4 |  -19.184  (>$10^3$)|   -5.006   (0.729) |
  TW:1          |  19.243   (>$10^3$)|    2.503   (0.532) |
  TW:2          |   2.111    (0.601) |    2.503   (0.532) |
  TW:3          |   2.940    (0.828) |    2.503   (0.532) |
  TW:4          |  17.064   (>$10^3$)|    2.503   (0.532) |
  CY:1          |   1.659    (1.179) |    1.528   (0.474) |  
  CY:2          |   1.343    (0.583) |    1.528   (0.474) |
  CY:3          |   1.693    (0.660) |    1.528   (0.474) |
  CY:4          |   1.162    (0.905) |    1.528   (0.474) |

Table: Estimates and standard errors (in parenthesis) of regression coefficients of the non-proportional odds model (NPOM) of the wine dataset, having temperature:warm (TW) and contact:yes (CY) as predictors, with estimates obtained using the `serp` and `vglm` R functions respectively penalized and non-penalized. 


As observed, with extreme shrinkage all subject-specific estimates (including the unbounded estimates) get smoothed out to parallel estimates (see the dashed blue line in each panel). Thus, estimates between the parameter space of NPOM and POM inclusive could be determined via an appropriate tuning procedure. In this instance, a 5-fold cross-validated tuning with SERP resulted in POM (see Table 1) with no boundary estimates. The intercept parameters too, with no penalty on them, were all stabilized by the smoothing penalty on the remaining coefficients. This obviously, demonstrates the amount of flexibility that could be achieved with `serp` in fitting cumulative models in general. In addition to the ordinal logit model, `serp` also fits ordinal models with the probit, loglog, cloglog and cauchit link functions. A penalized partial proportional odds model [@peterson_partial_1990] with the different link functions is also possible with `serp`.



# Conclusion
The R add-on package `serp` contains functions for regularization/smoothing across response categories in the non-proportional cumulative ordinal regression model. Beyond the highlighted functionality, `serp` provides a collection of tools that promote stress-free modeling in empirical research. Moreover, standard function names and arguments already known to users familiar with related libraries are also used in `serp`, reducing unnecessary ambiguity. Lastly, details about usage and more elaborate examples are hosted online through a pkgdown [@wickham_pkgdown_2020] website on [Github Pages](https://ejikeugba.github.io/serp).  


# Acknowledgements
The author would like to thank Jan Gertheiss for his invaluable suggestions and contributions towards the development of `serp`.


# References

