---
title: 'SERP: An R package for smoothing in ordinal regression'
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
affiliation: 1
orcid: 0000-0003-2572-0023
affiliations:
  - name: Department of Mathematics and Statistics, School of Economics and Social Sciences, Helmut Schmidt University, Hamburg, Germany
index: 1
date: 25 July 2021
bibliography: paper.bib
---  
  


# Summary 

The use of specialized methods for the analysis of categorical data has drastically increased in recent times [@agresti_categorical_2002]. Scientists very frequently use one form of the categorical response model or the other to study relationships between different variables of interest and also make predictions and/or inferences. For instance, several findings arising from clinical trials, social surveys, market research, etc, very much rely on results obtained from the application of models like: the generalized linear model [@nelder_generalized_1972; @McCullagh&Nelder1989], the cumulative link models (CLMs) [@mcCullagh_regression_1980] and similar models [see, e.g., @fahrmeir_multivariate_2001; @tutz_regression_2011]. However, in high-dimensional applications of these models where several potential predictors are used in the model, the so-called identification problem (particularly in CLMs) could be somewhat unavoidable [@fisher_identification_1966; @greenberg_advanced_1983]. Thankfully, the regularization technique [see, e.g., Hastie, Tibshirani, and Friedman 2009; @buhmann_statistics_2011] among other approaches, offers a means to address such an issue. 

The \texttt{serp} library [@Ugba_serp_2021], available in the Comprehensive R Archive Network ([CRAN](https://CRAN.R-project.org/package=serp)) [@R_language_2020], implements a unique regularization algorithm that results in very flexible modelling of CLMs [see, @tutz_regularized_2016; @ugba_smoothing_2021]. Fitting is powered by a modified-Newton optimization method that facilitates an easy convergence and estimation speed. As an open-source library, \texttt{serp} aims to provide a platform for advanced scientific modelling. The core features of \texttt{serp} as well as details on its usage are provided in the library [documentation](https://cran.r-project.org/web/packages/serp/serp.pdf).



# Statement of Need

Existing software implementations of some of the available regularization penalties in the literature, including: lasso [@tibshirani_regression_1996], ridge [@hoerl_ridge_1970] and elastic net [@zou_regularization_2005], seem to concentrate more on the linear and the binary outcome models than the ordinal models. For instance, the two popular CRAN packages for regularized regression, namely, the \texttt{glmnet} [@friedman_regularization_2010] and \texttt{penalized} [@goeman_penalized_2014] do not currently extend to ordinal models. Attempts to bridge this gap are, for instance, found in libraries such as \texttt{rms} [@harrell_regressionl_2015], \texttt{glmnetcr} [@archer_glmnetcr_2014a],  \texttt{glmpathcr} [@archer_glmpathcr_2014b], and \texttt{ordinalNet} [@wurm_regularized_2017], which all provide some means of regularizing the different forms of the ordinal model, using the aforementioned penalties. Similarly, the \texttt{serp} library, provides flexible modelling of CLMs that uses the smooth-effects-on-response penalty (SERP) [@ugba_smoothing_2021; @Ugba_serp_2021; @tutz_regularized_2016], a unique penalty that enforces smoothness on the adjacent categories of CLMs. In other words, when a pure non-proportional odds model (NPOM) is optimized together with SERP, one obtains the proportional odds model (POM) under extreme shrinkage. Unlike other commonly implemented regularization techniques which primarily act on the predictors (horizontal penalization), SERP provides a vertical penalization of the response categories, thus, bridging the gap between the general and the restricted form of CLMs. \texttt{serp} comes with several useful functions that support empirical research, being also equipped with standard model performance and descriptive methods.


# A minimal Example

The [wine dataset](https://ejikeugba.github.io/serp/reference/wine.html) adapted from @Randall_1989 represents the outcome of a factorial experiment on factors determining the bitterness of wine. Two treatment factors (temperature and contact) with two levels each (yes, no) are provided, with the rating of wine bitterness ranging from 1 = 'least bitter' to 5 = 'most bitter'. A typical modelling framework in this instance is an ordinal regression where, for instance, the cumulative logit model is used to predict the scores of wine bitterness using the two stated treatment factors. However, the fitted NPOM (see, \@ref(tab:estimates)) using, for instance, the vglm function from the R \texttt{VGAM} library [@yee_vgam_2010] is not fully identifiable. As observed, unbounded estimates with large standard errors ($SE > 10^3$) were obtained for the first and last categories of temperature, and also the last threshold (Intercept:4). In an actual sense, the absolute values of the estimates and the estimated standard errors for the unbounded parameters seems to diverge to $\infty$ as the stopping criteria of the iterative fitting procedure used become stricter and the number of iterations allowed increases. Existing software implementations of CLMs (without regularization) has no means of handling such issues but would rather suggest using a different modelling approach altogether. However, by incorporating SERP into the modelling framework, a fully identified model with bounded estimates could be achieved. This demonstrates the amount of flexibility achieved with \texttt{serp} in fitting cumulative models. The coefficient paths for the estimated model via \texttt{serp}, and for increasing values of the tuning parameter $[\lambda: 10^-3, 10^5]$, are shown in \autoref{fig:serpfig}.


![Estimated coefficient for the ordinal model of the wine data when using SERP. The thick lines on the top displays are the category-specific coefficients associated with the two predictors, under increasing values of $\lambda$ on the range ($10^-3, 10^5$). The dashed horizontal (blue) and vertical (red) lines denote the parallel estimates and the selected estimates based on the minimum deviance, respectively. The down displays further illustrate SERP's smoothing steps from the category-specific to the parallel estimates, the solid black, gray and dashed blue line strokes are NPOM, SERP and POM estimates, respectively; with the dashed black lines indicating SERP estimates chosen via minimum deviance.\label{fig:serpfig}](serp_fig.png){ width=70%}


Table: (\#tab:estimates) Estimates and standard errors (in parenthesis) of regression coefficients of the non-penalized (via vglm) and penalized (via serp) NPOM of the wine dataset, with Temperature-warm (TW) and Contact-yes (CY) as predictors.


  Coefficients  |        vglm        |        serp        |
  :-------------|-------------------:|-------------------:|
  (Intercept):1 |   1.226    (0.557) |    1.255   (0.556) |
  (Intercept):2 |  -1.033    (0.481) |   -1.061   (0.480) | 
  (Intercept):3 |  -3.946    (0.902) |   -3.926   (0.881) |
  (Intercept):4 | -19.184   (>$10^3$)|   -6.497   (2.486) |
  TW:1          |  19.243   (>$10^3$)|    3.922   (2.333) |
  TW:2          |   2.111    (0.601) |    2.149   (0.600) |
  TW:3          |   2.940    (0.828) |    2.927   (0.810) |
  TW:4          |  17.064   (>$10^3$)|    4.368   (2.415) |
  CY:1          |   1.659    (1.179) |    1.624   (1.126) |  
  CY:2          |   1.343    (0.583) |    1.367   (0.580) |
  CY:3          |   1.693    (0.660) |    1.689   (0.652) |
  CY:4          |   1.162    (0.905) |    1.158   (0.892) |



Finally, beyond the discussed regularization functionality, \texttt{serp} provides a collection of tools that promote a stress-free modelling in empirical research. Moreover, standard function names and arguments already known to users familiar with related libraries are also used in \texttt{serp}, reducing unnecessary ambiguity. Lastly, details about usage and a more elaborate examples are hosted online through a pkgdown [@wickham_pkgdown_2020] website on [Github Pages](https://ejikeugba.github.io/serp).  


# Acknowledgements
A special thanks to Prof Dr. Jan Gertheiss for his invaluable suggestions and contributions towards the development of \texttt{serp}.


# References
