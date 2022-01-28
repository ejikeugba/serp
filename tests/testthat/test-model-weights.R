library(testthat)
library(serp)

context("check model weights and the errorMetrics function")
wine <- serp::wine

## checks on weight arguments
test_that("weight argument introduces no error",
          {
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "parallel", reverse=FALSE, weights = c(rep(1,71), NA),
                   data = wine))
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "parallel", reverse=FALSE, weights =  c(rep(1,71), -1),
                   data = wine))
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "parallel", reverse=FALSE, weights = c(rep(1,71), 0.6),
                   weight.type = "frequency", data = wine))
            expect_vector(
              serp(rating ~ temp + contact, slope = "parallel",
                   link = "cloglog", reverse=TRUE, weights = rep(1, 72),
                   weight.type = "frequency", data = wine)$coef)
            expect_vector(
              serp(rating ~ temp + contact, slope = "parallel",
                   link = "cloglog", reverse=TRUE, weights = rep(0.1, 72),
                   weight.type = "analytic", data = wine)$coef)

            set.seed(1)
            n <- 30
            test_data <- data.frame(y= as.ordered(rbinom(n,5, 0.1)),
                                    x1=runif(n), x2=rexp(n))
            expect_error(
              serp(y ~ x1 + x2 ,
                   slope = "penalize",
                   link = "logit",
                   tuneMethod = "cv",
                   gridType = "discrete",
                   weights = runif(nrow(test_data)),
                   reverse = F,
                   data = test_data)
            )

            rm(test_data, n)
          })

## checks on trace and errorMeterics
test_that("trace and errorMeterics work properly",
          {
            ## check if trace works
            expect_output(serp(rating ~ temp + contact, link = "logit",
                               slope = "unparallel", reverse=FALSE,
                               control= list(trace=1),
                               data=wine, subset = c(1:30)))

            expect_output(serp(rating ~ temp + contact, link = "loglog",
                               slope = "penalize", reverse=TRUE,
                               gridType = "fine",
                               control= list(trace=2),
                               data=wine, subset = c(1:50)))
            expect_output(serp(rating ~ temp + contact, link = "cloglog",
                               slope = "partial", reverse=TRUE,
                               globalEff = ~ temp + contact,
                               control= list(trace=3),
                               data=wine, subset = c(1:30)))

            ## checks on errorMetrics
            f1 <- serp(rating ~ temp + contact, link = "logit",
                       slope = "parallel", reverse=FALSE,
                       data = wine)
            hh <- list()
            hh$minp <- 1e-02
            fv <- f1$fitted.values
            expect_error(
              serp:::errorMetrics(f1$model[,1L], fv[,-1L], control = hh,
                                  type = "brier"))
            expect_vector(serp:::errorMetrics(f1$model[,1L], fv, control = hh, type = "logloss"))
            expect_vector(serp:::errorMetrics(f1$model[,1L], fv, control = hh, type = "misclass"))

            set.seed(1)
            y <- sample(c(0,1), 50, replace = TRUE)
            mm <- glm(y ~ rnorm(50))
             expect_error(serp:::vcov.serp(mm))
})
