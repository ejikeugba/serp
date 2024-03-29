
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable state and is being activelydeveloped](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test coverage](https://codecov.io/gh/ejikeugba/serp/branch/master/graph/badge.svg)](https://codecov.io/gh/ejikeugba/serp?branch=master)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/serp)](https://CRAN.R-project.org/package=serp)
[![CRAN status](https://www.r-pkg.org/badges/version/serp )](https://CRAN.R-project.org/package=serp)
[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0.en.html)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/ejikeugba/serp?branch=master&svg=true)](https://ci.appveyor.com/project/ejikeugba/serp)
[![R build status](https://github.com/ejikeugba/serp/workflows/R-CMD-check/badge.svg)](https://github.com/ejikeugba/serp/actions)
[![status](https://joss.theoj.org/papers/6ebd3b75ea792be908f0dadebd7cf81c/status.svg)](https://joss.theoj.org/papers/6ebd3b75ea792be908f0dadebd7cf81c)
<!-- badges: end -->

# Contributing to serp

Thank you for taking the time to contribute to the development of `serp`. You could find the following guidelines useful in making your contributions. 

Before you start:

* It is important to have a valid [GitHub account](https://github.com/signup/free).
* Trivial changes to comments or documentation do not require creating a new issue.

## Did you find a bug?

* Make sure the bug was not already reported in the Github [Issues](https://github.com/ejikeugba/serp/issues).
* [Open an issue](https://github.com/ejikeugba/serp/issues/new) and clearly describe the issue with as much information as possible. A code sample or an executable test case are recommended.
  
## Did you plan to write a patch that fixes a bug?

  * [Open an issue](https://github.com/ejikeugba/serp/issues/new) and clearly describes the problem and discuss how your solution will affect `serp`.
  * Fork the repository on GitHub to work on the patch.
  * Get in touch with the maintainer to refine and prioritize your issue.

## Making changes and Pull requests

* Start your work on your fork of the repository. If you haven't done this before, try using `usethis::create_from_github("ejikeugba/serp", fork = TRUE)`.
* Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
* Create a Git branch for your pull request (PR). You may want to use `usethis::pr_init("brief-description-of-change")`.
* Check for unnecessary whitespace with `git diff --check` and format code.
* Commit messages should be descriptive, mentioning what was changed and why, and also **reference the relevant issue number**. 
* Ensure to add the necessary tests for your changes (testthat preferably).
* Run **all** the tests to assure nothing else was accidentally broken, also keep an eye on the test coverage metric.
* Commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser
    
## Copyright issues
 
* On submission, it is crucial your PR includes the following statement: You own the copyright on the code being contributed, and you hereby grant `serp` repo cph unlimited license to use this code in this version or any future version of `serp`. You reserve all other rights to the code.
* It may not be advisable to contribute third party codes to this project. Useful suggestions are nonetheless welcomed.
* The Pull Requests are thereafter reviewed, with feedbacks communicated as soon as possible.

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [About pull requests](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)
