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
affiliation: 1
orcid: 0000-0003-2572-0023
affiliations:
  - name: Department of Mathematics and Statistics, School of Economics and Social Sciences, Helmut Schmidt University, Hamburg, Germany
index: 1
date: 25 July 2021
bibliography: paper.bib
---  
  


# Summary 

The use of specialized methods for the analysis of categorical data has drastically increased in recent times [@agresti_categorical_2002]. Scientists frequently use the different forms of the ordinal model to define relationships between study variables of interest and for prediction/inference. For instance, the cumulative link model (CLM) propagated by @mcCullagh_regression_1980 finds a wide range of empirical applications in clinical trials, social surveys, market research, etc. However, in high-dimensional settings where several potential predictors and response categories are modelled, the so-called identification problem could be somewhat unavoidable [@fisher_identification_1966; @greenberg_advanced_1983]. Thankfully, the regularization technique [see, e.g., Hastie, Tibshirani, and Friedman 2009; @buhmann_statistics_2011] among other approaches, offers a remedy in such situations.    

\texttt{serp} software package [@Ugba_serp_2021], available in the Comprehensive R Archive Network ([CRAN](https://CRAN.R-project.org/package=serp)) [@R_language_2021], implements a unique regularization algorithm that enforces smoothness on the adjacent categories of CLMs. The approach provides flexible modeling of CLMs that uses the smooth-effects-on-response penalty (SERP) discussed in @ugba_smoothing_2021, [see also, @tutz_regularized_2016]. When applied to the general form of the cumulative logit model, otherwise known as the non-proportional odds model (NPOM), one obtains the proportional odds model (POM) under extreme parameter shrinkage. Unlike most implemented regularization penalties for CLMs that primarily act on the predictors, a form of horizontal penalization, SERP provides a vertically oriented penalization that focuses on the response categories, thus, making it possible to regularize the parameter space between the general and restricted form of the cumulative link models. Fitting in \texttt{serp} is powered by a modified Newton's optimization method that facilitates an easy convergence and estimation speed. Moreover, as an open-source software package, \texttt{serp} aims to provide a platform for advanced scientific modeling in empirical research. The core features of \texttt{serp}, as well as details about its usage, are provided in the software package [documentation](https://cran.r-project.org/web/packages/serp/serp.pdf).



# Statement of Need

Although a very useful tool in regression analysis, the development of regularization methods for the ordinal categorical response models is still in its elementary stage. This explains why so many popular statistical software packages for ordinal regression do not currently implement such methods. For instance, the vglm function in the \texttt{VGAM} R package [@yee_vector_1996; @yee_vgam_2010] fits the cumulative family of the ordered model, with several link functions, but with no provision for regularization. The same is true with the SAS \texttt{CATMOD} procedure [@SAS_Institute_2018]. Exiting software packages so far implementing regularization techniques for ordinal models are mainly extensions of methods typically found in the high-dimensional linear or generalized linear model framework, where penalties in the tradition of lasso [@tibshirani_regression_1996], ridge [@hoerl_ridge_1970] and elastic net [@zou_regularization_2005] are mostly implemented. For instance, the glmnetcr [@archer_glmnetcr_2014a] and glmpathcr [@archer_glmpathcr_2014b] R packages, fit continuation ratio models with the elastic net penalty, while the rms R package [@harrell_regressionl_2015] fits the cumulative logit model with quadratic (ridge regression) penalty. Also, the \texttt{ordinalgmifs} R package provided by @archer_ordinalgmifs_2014 implements a generalized monotone incremental forward stagewise (GMIFS) algorithm for regularized ordinal regression models, with a solution path akin to the L1 norm (lasso) penalty. Furthermore, the \texttt{ordinalNet} R package [@wurm_regularized_2017] fits the elementwise link multinomial-ordinal (ELMO) class of models with the elastic net using a coordinate descent algorithm.

A common principle behind the regularization procedures with the highlighted software packages is to make the model's maximum likelihood estimates a bit biased towards zero or fully shrunk to zero, of course, with the so-called bias-variance trade-off in mind. However, the regularization method implemented in \texttt{serp} provides a complete paradigm shift from business as usual. \texttt{serp}  provides a unique form of penalization that shrinks the subject-specific effects in the cumulative link model towards global effects under extreme shrinkage. This form of penalization results in the proportional odds model in terms of the cumulative logit model. A somewhat close form of penalization that is based on smoothing estimates in categorical response models is provided in the \texttt{mgcv} R package [@wood_smoothing_2016] but for un-ordered categories and with smoothing done in the linear predictors. By and large, \texttt{serp} provides a regularization routine for smoothing parameter estimates in the cumulative ordinal models, with resultant reduction in model complexity that doesn't require chopping off important features in the model. The \texttt{serp} software package comes with several useful functions that support empirical research, being also equipped with standard model performance and descriptive methods.



# A minimal Example

The [wine dataset](https://ejikeugba.github.io/serp/reference/wine.html) adapted from @Randall_1989 represents the outcome of a factorial experiment on factors determining the bitterness of wine. Two treatment factors (temperature and contact) with two levels each (yes, no) are provided, with the rating of wine bitterness ranging from 1 = 'least bitter' to 5 = 'most bitter'. A typical modeling framework in this instance is an ordinal regression where, for instance, the cumulative logit model is used to predict the scores of wine bitterness using the two stated treatment factors. However, the fitted NPOM (see, Table 1) using the \texttt{VGAM}-vglm function [@yee_vgam_2010], for instance,  is not fully identifiable. As observed, unbounded estimates with large standard errors ($SE > 10^3$) were obtained for the first and last categories of temperature, and also the last threshold (Intercept:4). In a general sense, the absolute values of the estimates and the estimated standard errors for the unbounded parameters seem to diverge to $\infty$ as the stopping criteria of the iterative fitting procedure used become stricter and the number of iterations allowed increases. Other software implementations of the same model without regularization run into a similar problem, forcing users to adopt a different modeling approach which may or may not be appropriate for the data. However, by incorporating SERP into the modeling framework, a fully identified model with bounded estimates could be achieved (as shown in Table 1). This obviously, demonstrates the amount of flexibility that could be achieved with the \texttt{serp} software in fitting the cumulative models in general. The coefficient paths for the estimated model via \texttt{serp}, using increasing values of the tuning parameter $[\lambda: 10^-3, 10^5]$, are shown in \autoref{fig:serpfig}.


![Estimated coefficient paths for the ordinal model of the wine data when using SERP. The thick lines on the top displays are the category-specific coefficients associated with the two predictors, under increasing values of $\lambda$ on the range ($10^-3, 10^5$). The dashed horizontal (blue) and vertical (red) lines denote the parallel estimates and the selected estimates based on the minimum deviance, respectively. The down displays further illustrate SERP's smoothing steps from the category-specific to the parallel estimates, the solid black, grey and dashed blue line strokes are NPOM, SERP and POM estimates, respectively; with the dashed black lines indicating SERP estimates chosen via minimum deviance.\label{fig:serpfig}](serp_fig.png){ width=70%}



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


Table: Estimates and standard errors (in parenthesis) of regression coefficients of the non-penalized (via vglm) and penalized (via serp) NPOM of the wine dataset, with Temperature:warm (TW) and Contact:yes (CY) as predictors.


# Conclusion
Beyond the highlighted regularization functionality, \texttt{serp} provides a collection of tools that promote stress-free ordinal regression modeling in empirical research. Moreover, standard function names and arguments already known to users familiar with related libraries are also used in \texttt{serp}, reducing unnecessary ambiguity. Lastly, details about usage and more elaborate examples are hosted online through a pkgdown [@wickham_pkgdown_2020] website on [Github Pages](https://ejikeugba.github.io/serp).  


# Acknowledgements
A special thanks to Prof Dr. Jan Gertheiss for his invaluable suggestions and contributions towards the development of \texttt{serp}.


# References

