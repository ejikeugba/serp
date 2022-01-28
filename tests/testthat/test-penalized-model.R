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
            #1# From NPOM to POM
            sp1 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "penalize", tuneMethod = "user",
                        lambda = 1e07, reverse=FALSE, data=wine,
                        subset = c(1:30))
            cof1 <- round(as.numeric(coef(sp1))[1L], 4L)

            sp2 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "parallel", reverse=FALSE, data = wine,
                        subset = c(1:30))
            cof2 <- round(as.numeric(coef(sp2))[1L], 4L)
            expect_equal(cof1, cof2, check.attributes=FALSE,
                         tolerance=tol)

            sp3 <- serp(rating ~ temp + contact, link = "probit",
                        slope = "penalize", tuneMethod = "cv",
                        lambdaGrid = 10^seq(-1, 2, length.out=2),
                        gridType = "fine", reverse=TRUE, data=wine,
                        control = list(nrFold =5),
                        subset = c(1:30))

            sp4 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "penalize", tuneMethod = "finite",
                        lambda = 1e10, reverse=FALSE, data=wine,
                        subset = c(1:50))

            expect_output(penalty.print(object=sp1, max.tun=TRUE))
            expect_output(penalty.print(object=sp3, max.tun=TRUE))
            expect_output(penalty.print(object=sp4, max.tun=TRUE))
            expect_vector(summary(sp3)$penalty$lambda)
            expect_null(summary(sp2)$penalty)

            #2# partial slope with all variables in globalEff yields  POM
            sp5 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "partial", reverse=FALSE, subset = c(1:50),
                        globalEff = ~ temp + contact, data=wine)
            sp6 <- serp(rating ~ temp + contact, link = "logit",
                        slope = "parallel", reverse=FALSE, subset = c(1:50),
                        data=wine)

            #3# checks on anova, confinct and vcov functions
            expect_equal(coef(sp5), coef(sp6), check.attributes=FALSE,
                         tolerance=tol)
            expect_error(anova.serp(sp5))
            expect_false(inherits(try(anova.serp(sp5, sp6)), 'try-error'))
            expect_error(anova.serp(lm(rnorm(50) ~ runif(50))))
            expect_error(anova.serp(sp3, sp4))
            expect_error(anova.serp(sp5, update(sp6, subset = c(1:40))))
            expect_vector(length(anova.serp(sp5, sp6, test = "none")))
            expect_error(confint.serp(sp5, sp6))
            expect_false(inherits(try(confint.serp(sp5)), 'try-error'))
            expect_error(confint.serp(sp5, level = 1.2))
            expect_message(confint.serp(sp5, parm = 0.1))
            expect_error(confint.serp(lm(rnorm(50) ~ runif(50))))
            expect_false(inherits(try(vcov.serp(sp5)), 'try-error'))
            expect_error(vcov.serp(sp1, sp2))
            expect_vector(AIC.serp(sp5))
            expect_vector(BIC.serp(sp5))
            expect_vector(logLik.serp(sp5))
            expect_output(print.serp(sp5))
            expect_output(print.summary.serp(summary.serp(sp5)))

            rm(sp1, sp2, sp3, sp4, sp5, sp6, cof1, cof2)
          })

## checks on lambda and lambdaGrid
test_that("lambda is a single numeric and non-negative value and
           that lambdaGrid contains the right input values",
          {
            subs <- c(1:30)
            expect_error(
              serp(rating ~ temp + contact, slope = "penalize",
                   link = "loglog", reverse = TRUE, tuneMethod = "user",
                   lambda = c(0.3,5), data = wine, subset = subs))

            expect_error(
              serp(rating ~ temp + contact, link = "probit",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "aic",
                   lambdaGrid = c(-3,4),
                   data = wine, subset = subs))

            expect_vector(
              serp(rating ~ temp + contact, link = "cloglog",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "aic",
                   lambdaGrid = c(0,10),
                   data = wine, subset = subs)$lambda)

            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse = FALSE,
                   tuneMethod = "deviance",
                   lambdaGrid = c(0,1), control = list(trace=-1),
                   data = wine, subset = subs))
            expect_error(serp.control(misclass.thresh = 5))
            expect_error(serp.control(maxpen = 1e11))
            rm(subs)
          })

## checks on data subset
test_that("subset indices are positive whole numbers",
          {
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, "r"),
                   data = wine))
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, 4.7),
                   data = wine))
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, -5),
                   data = wine))
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "penalize", reverse=TRUE, subset = c(1, 2, 3, NA),
                   data = wine))
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

            s1 <- serp(rating ~ temp + contact, link = "logit",
                       slope = "penalize", reverse = TRUE, gridType = "fine",
                       globalEff = ~ contact, data = wine, subset = subs)
            pred2 <- predict(s1, type = 'link', newdata=head(wine))
            expect_vector(pred2[1L,])

            wdat <- head(wine)
            colnames(wdat) <- 1:6

            expect_error(predict(s1, type = 'link',
                                 newdata = wdat))

            rtt <-  wine$rating
            rtt <- data.frame(rating=rtt)

            s2 <- serp(rating ~ 1, link = "logit", slope = "penalize",
                       reverse = TRUE, gridType = "fine",
                       data =rtt)
            pred3 <- predict(s2, type = 'link', newdata=head(wine))
            expect_vector(pred3[1L,])

            rm(pred1, pred2, pred3, subs, wdat, rtt, s2)

          })

## checks on warning and error messages
test_that("error messages and warnings report properly",
          {set.seed(1)
            n <- 8
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
                   data = test_data1))

            expect_error(print.serp(test_data1))
            expect_error(print.summary.serp(test_data1))
            expect_error(summary.serp(test_data1))
            expect_error(predict.serp(test_data1))

            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "cv",
                   control = list(nrFold = 1),
                   data = test_data1))


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
                   data = test_data2))


            set.seed(11)
            n <- 20
            test_data3 <- data.frame(y= as.ordered(rbinom(n,3, 0.2)),
                                     x1=runif(n), x2=runif(n))
            expect_warning(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "finite",
                   gridType = "fine",
                   reverse = F,
                   lambdaGrid = 10^seq(-1, 2, length.out=2),
                   data = test_data3))


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
                           lambdaGrid = 10^seq(-1, 2, length.out=2),
                           control = list(nrFold =2),
                           data = test_data4))

            expect_output(penalty.print(object=mm4, max.tun=TRUE))
            expect_output(print.serp(mm4))

            expect_error(
              serp(y ~ x1 + x2,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "user",
                   reverse = TRUE,
                   lambda,
                   data = test_data4))

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
                   data = test_data4))

            expect_error(
              serp(rating ~ temp * contact,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "user",
                   data = wine))

            expect_error(
              serp(rating ~ temp * contact,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "user",
                   lambda = 1e11,
                   data = wine))

            expect_error(
              serp(rating ~ temp * contact,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "aic",
                   lambdaGrid = 10^seq(-1, 12, length.out=2),
                   data = wine))

            expect_error(
              serp(rating ~ temp * contact, slope = "partial",
                   link = "cauchit", globalEff= ~1, data = wine))

            sdat <- wine
            sdat$extra <- 1:72
            expect_error(
              serp(rating ~ temp * contact, slope = "partial",
                   link = "cauchit", globalEff= ~extra, data = sdat))

            test_data4$y <- rbinom(20, 1, 0.5)
            expect_error(
              serp(ordered(y) ~ x1 + x2,
                   slope = "parallel",
                   link = "logit", data = test_data4))

            rm(test_data1, test_data2, test_data3, test_data4, n)
          })
