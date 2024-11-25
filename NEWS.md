## serp 0.2.5
- Resolved issues causing failures in CRAN checks on multiple platforms.
- Added `tibble`, `vctrs`, and `pkgdown` to `Suggests`.
- Minor fixes applied to the `serp.fit` file.

---

## serp 0.2.4
- The `errorMetrics()` function is no longer exported. An equivalent function with enhanced functionality is provided by the `erroR()` function in the **gofcat** package.
- Fixed dependency-induced bugs in package tests.
- Removed the `License` file from the `DESCRIPTION` and associated file docs.
- Minor updates to the **serp** documentation.

---

## serp 0.2.3
- CRAN release.
- Imported `crayon` to enable colored outputs in returned objects.
- Improved test coverage by updating **serp** tests.
- Fixed a bug in `errorMetrics`, with the `model` argument removed.
- Added a `print` method and object class for `errorMetrics` results.
- Replaced `cat()` with `message()` where appropriate.
- Added a citation for the **serp** package.

---

## serp 0.2.2
- JOSS release.
- Minor updates to **serp** documentation.
- Provided `coefficient()` as an alias for `coef()`.
- Updated `README.md` to include community guidelines and a contributors' code of conduct.

---

## serp 0.2.1
- Initial submission to CRAN.

---

## serp 0.2.0.9001
- Added the `update` function to the namespace.
- Included examples in the documentation of various functions.
- Set the shrinkage parameter's upper limit to `1e10` in `serp.control`.
- Fixed bugs in the test functions.

---

## serp 0.2.0.9000
- Replaced the deviance tuning option in `serp.tuneMethod` with AIC.
- **serp** output now includes residual degrees of freedom (rdf).
- Added functionality to compute effective degrees of freedom and rdf from the trace of the generalized hat matrix.
- Updated the "Value" section of the `serp.summary` documentation.
- Enhanced **serp** test functions.
- Made changes to `README.Rmd` and `README.md`.

---

## serp 0.2.0
- Released version 0.2.0 of **serp**.

---

## serp 0.1.9.9001
- Re-submitted to CRAN.

---

## serp 0.1.9.9000
- Fixed all issues identified in the last CRAN release.

---

## serp 0.1.9
- Submitted **serp** version 0.1.9 to CRAN.

---

## serp 0.1.8.9007
- Updated `README.md`.
- Added DOI to the `DESCRIPTION` file.
- Updated references in **serp** documentation.

---

## serp 0.1.8.9006
- Fixed bugs in tests.

---

## serp 0.1.8.9005
- Reconstructed tests to improve unit test coverage.
- Revised error and warning messages in key **serp** functions.
- Fixed bugs in `anova.serp`.

---

## serp 0.1.8.9004
- Updated the `License` field in `DESCRIPTION`.
- Fixed bugs in example code.
- Corrected bugs in test files.

---

## serp 0.1.8.9003
- Corrected column titles for predicted values with `reverse = TRUE`.
- Fixed bugs in the deviance function (`dvfun`).
- Enhanced the `serp.predict` function to handle single and multiple row inputs.
- Improved the iteration algorithm in `serp.fit`.
- Removed commented lines from the main **serp** function.
- Moved `reverse` and `linkf` arguments from the main function to `serp.fit`.
- Corrected results for `loglog` and `cloglog` links with `reverse = TRUE`.
- Split long lines of code in the score function for readability.
- Updated starting values to account for different link functions.
- Removed the `reverse` statement from the `cv` function.
- Deprecated `reverse.fun` in `summary.serp`.
- Normalized predicted values in `errorMetrics` for efficiency.

---

## serp 0.1.8.9002
- Updated warning messages in `serp.fit`.
- Fixed a bug in `serp.fit` when specifying `gridType`.
- Removed `trainError` from **serp** output.
- Deprecated `reverse.fun`.
- Made minor adjustments to `penalty.print`.
- Fixed bugs in `print.summary.serp`.
- Adjusted `summary.serp` to exclude `trainError`.
- Reconstructed the `predict` function and fixed bugs related to the `newdata` argument.
- Fixed bugs in the `anova` function, now including slope type in the output.
- Disabled `anova` tests for the `penalize` slope type.
- Updated warning messages in `serp.control`.

---

## serp 0.1.8.9001
- Adopted a new license for **serp**.

---

## serp 0.1.8.9000
- Added badges and a hexagon sticker to `README.md`.

---

## serp 0.1.8
- Re-submitted the **serp** package to CRAN.
- Included references to methods used in the package in the `DESCRIPTION` field.
- Added `Value` sections to all exported methods with detailed explanations.
- Enabled all examples by removing `\dontrun{}` blocks.
- Removed all occurrences of `<<-` in functions.
- Updated `README.md`.
