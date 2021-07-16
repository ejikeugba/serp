library(testthat)
library(serp)

context("Penalized - checks if serp works properly on cumulative models")
wine <- serp::wine

test_that("estreem shrinkage with serp results to the
           proportional odds model (POM), while zero
           shrinkage returns the non-proportional odds
           model (NPOM)",
          {
            tol <- 1e-06
            subs <- c(1:30)
            #1# From NPOM to POM
            sp1 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "penalize", tuneMethod = "user",
                        lambda = 1e10, reverse=FALSE, data=wine,
                        subset = subs)
            cof1 <- unique(round(as.numeric(coef(sp1)), 9L))

            sp2 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "parallel", reverse=FALSE, data = wine,
                        subset = subs)
            cof2 <- coef(sp2)
            expect_equal(cof1, cof2, check.attributes=FALSE,
                         tolerance=tol)

            sp3 <- serp(rating ~ temp + contact, link = "probit",
                        slope = "penalize", tuneMethod = "cv",
                        lambdaGrid = 10^seq(-1, 2, length.out=2),
                        gridType = "fine", reverse=TRUE, data=wine,
                        subset = subs)

            sp4 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "penalize", tuneMethod = "finite",
                        lambda = 1e10, reverse=FALSE, data=wine,
                        subset = c(1:50))

            expect_output(penalty.print(object=sp1, max.tun=TRUE))
            expect_output(penalty.print(object=sp3, max.tun=TRUE))
            expect_output(penalty.print(object=sp4, max.tun=TRUE))
            expect_vector(summary(sp3)$penalty$lambda)
            expect_null(summary(sp2)$penalty)

            rm(sp1, sp2, cof1, cof2, subs)

            #2# partial slope with all variables in glodalEff yields  POM
            sp1 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "partial", reverse=FALSE, subset = c(1:50),
                        globalEff = ~ temp + contact, data=wine)
            sp2 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "parallel", reverse=FALSE, subset = c(1:50),
                        data=wine)

            #3# checks on anova, confinct and vcov functions

            expect_equal(coef(sp1), coef(sp2), check.attributes=FALSE,
                         tolerance=tol)
            expect_error(anova.serp(sp1),
                         "no anova implementation yet for a single 'serp' object")

            expect_false(inherits(try(anova.serp(sp1, sp2)), 'try-error'))

            expect_error(confint.serp(sp1, sp2), "one object at a time allowed")

            expect_false(inherits(try(confint.serp(sp1)), 'try-error'))

            expect_error(vcov.serp(sp1, sp2), "one object at a time allowed")

            expect_false(inherits(try(vcov.serp(sp1)), 'try-error'))

            expect_vector(AIC.serp(sp1))
            expect_vector(BIC.serp(sp1))
            expect_vector(logLik.serp(sp1))
            expect_output(print.serp(sp1))
            expect_output(print.summary.serp(summary.serp(sp1)))

            rm(sp1, sp2, sp3, sp4)
          })

## checks on lambda and lambdaGrid
test_that("lambda is a single numeric and non-negative value and
           that lambdaGrid contains the right input values",
          {
            subs <- c(1:30)
            expect_error(
              serp(rating ~ temp + contact, slope = "penalize",
                   link = "loglog", reverse = TRUE, tuneMethod = "user",
                   lambda = c(0.3,5), data = wine, subset = subs),
              "lambda should be a single numeric and non-negative value")

            expect_error(
              serp(rating ~ temp + contact, link = "probit",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "deviance",
                   lambdaGrid = c(-3,4),
                   data = wine, subset = subs),
              "lambdaGrid must be a non-negative numeric vector of length > 1")

            expect_vector(
              serp(rating ~ temp + contact, link = "cloglog",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "deviance",
                   lambdaGrid = c(0,10),
                   data = wine, subset = subs)$lambda)

            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "deviance",
                   lambdaGrid = c(0,1), control = list(trace=-1),
                   data = wine, subset = subs),
              "maxits, eps, maxpen, minP, maxAdjIter, max.half.iter and
  relTol should all be numeric and non-negative")
            rm(subs)
          })

## checks on data subset
test_that("subset indices are positive whole numbers",
          {
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, "r"),
                   data = wine),
              "subset indices must be positive whole numbers")
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, 4.7),
                   data = wine),
              "subset indices must be positive whole numbers")
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, -5),
                   data = wine),
              "subset indices must be positive whole numbers")
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, NA),
                   data = wine),
              "subset indices must be positive whole numbers")
          })

## checks on predict function and data subset
test_that("predict function works properly",
          {
            subs <- c(1:30)
            pred1 <- predict(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse = FALSE, gridType = "fine",
                   globalEff = ~ temp + contact, data = wine, subset = subs),
              type = 'link', newdata=head(wine))

            expect_vector(pred1[1L,])

            pred2 <- predict(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse = TRUE, gridType = "fine",
                   globalEff = ~ temp + contact, data = wine, subset = subs),
              type = 'link', newdata=head(wine))

            expect_vector(pred2[1L,])

            expect_error(
              errorMetrics(
                actual = wine$rating,
                predicted = matrix(runif(100), 20,5),
                model = c("multiclass"),
                type = c("brier"),
                eps = .Machine$double.eps),
              "levels of actual observations not equal to the number of columns of fitted values, or unequal lengths of observations")

            expect_error(
              errorMetrics(
                actual = as.factor(rbinom(100,1,0.5)),
                predicted = matrix(runif(100), 20,1),
                model = c("binary"),
                type = c("brier"),
                eps = .Machine$double.eps),
              "supply a vector of fitted values")

            expect_error(
              errorMetrics(
                actual = as.factor(rbinom(100,1,0.5)),
                predicted = c(runif(98),'x','y'),
                model = c("binary"),
                type = c("brier"),
                eps = .Machine$double.eps),
              "supply a numeric vector of fitted values")


            rm(pred1, pred2, subs)

          })

## checks on warning and error messages
test_that("error messages and warnings report properly",
          {set.seed(1)
            n <- 10
            test_data1 <- data.frame(y= as.ordered(rbinom(n,3, 0.2)),
                                     x1=rnorm(n), x2=runif(n))
            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "cv",
                   gridType = "discrete",
                   reverse = F,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data1),
              "cv tuning did not succeed, try switching to a different tuneMethod")

            expect_error(print.serp(test_data1),
                         "input must be an object of class 'serp'")
            expect_error(print.summary.serp(test_data1),
                         "input must be an object of class 'serp'")
            expect_error(summary.serp(test_data1),
                         "input must be an object of class 'serp'")
            expect_error(predict.serp(test_data1),
                         "object must be of class \"serp\"")

            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "cv",
                   control = list(nrFold = 1),
                   data = test_data1),
              "nrFold should be numeric and between 2 and 10 inclusive.")


            set.seed(1)
            n <- 20
            test_data2 <- data.frame(y= as.ordered(rbinom(n,3, 0.2)),
                                     x1=runif(n), x2=runif(n))
            expect_warning(
              serp(y ~ x1 + x2,
                   slope = "unparallel",
                   link = "logit",
                   tuneMethod = "cv",
                   gridType = "discrete",
                   reverse = TRUE,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data2),
              "Stochastic ordering assumption failed.
    Consider using the penalized, parallel or partial slope,
    or other link functions.")


            set.seed(1)
            n <- 20
            test_data3 <- data.frame(y= as.ordered(rbinom(n,3, 0.1)),
                                     x1=runif(n), x2=runif(n))
            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "cv",
                   gridType = "fine",
                   reverse = F,
                   lambda = 0.09,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data3),
              "cv tuning did not succeed, try switching to a different tuneMethod")

            set.seed(1)
            n <- 20
            test_data4 <- data.frame(y= as.ordered(rbinom(n,2, 0.3)),
                                     x1=runif(n), x2=rexp(n))

            expect_warning(
              mm4 <-  serp(y ~ x1 + x2,
                           slope = "penalize",
                           link = "logit",
                           tuneMethod = "cv",
                           gridType = "discrete",
                           reverse = F,
                           lambda = 0,
                           lambdaGrid = 10^seq(-1, 2, length.out=2),
                           data = test_data4),
              "non-finite log-likelihood persists, increase lambdaGrid upper limit or apply a different tuning method")

            expect_output(penalty.print(object=mm4, max.tun=TRUE))
            expect_output(print.serp(mm4))

            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "user",
                   gridType = "discrete",
                   reverse = TRUE,
                   lambda,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data4),
              "unassigned value in serp function.")

            expect_vector(
              serp(y ~ 1,
                   slope = "unparallel",
                   link = "logit",
                   tuneMethod = "user",
                   gridType = "discrete",
                   reverse = TRUE,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data4)$logLik)

            test_data4$y <- rpois(20, 1)
            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "user",
                   gridType = "discrete",
                   reverse = TRUE,
                   lambda = 4,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data4),
              "response must be an ordered factor")

            rm(test_data1, test_data2, test_data3, test_data4, n)
          })
