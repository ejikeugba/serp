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
                   data = wine),
              "weights should be numeric vector with no NA's")
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "parallel", reverse=FALSE, weights =  c(rep(1,71), -1),
                   data = wine),
              "negative weights not allowed")
            expect_error(
              serp(rating ~ temp + contact, link = "logit",
                   slope = "parallel", reverse=FALSE, weights = c(rep(1,71), 0.6),
                   weight.type = "frequency", data = wine),
              "frequency weights must be whole numbers")
            expect_vector(
              serp(rating ~ temp + contact, slope = "parallel",
                   link = "cloglog", reverse=TRUE, weights = rep(1, 72),
                   weight.type = "frequency", data = wine)$coef)

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
                   data = test_data),
              "only frequency weights are allowed in 'cv' tuning."
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
            expect_error(
              errorMetrics(f1$model[,1L], f1$model[,1L],type = "brier"),
              "supply either a matrix or dataframe of fitted values")

            e1 <- errorMetrics(f1, type = "brier")
            e2 <- errorMetrics(f1, type = "logloss")
            e3 <- errorMetrics(f1, type = "misclass")

            ## set.seed(1)
            y <- sample(c(0,1), 50, replace = TRUE)
            mm <- glm(y ~ rnorm(50))
            p <- mm$fitted.values

            expect_error(
              errorMetrics(mm$y, mm$fitted.values, model= "binary", type = "brier"),
              "'actual' must be a factor")

            yna <- factor(c(y,1))
            pna <- c(p, NA)
            y <- factor(mm$y)

            e4 <- errorMetrics(y, p, model= "binary", type = "logloss")
            e5 <- errorMetrics(y, p, model= "binary", type = "misclass")

            expect_false(any(is.na(c(e1, e2, e3, e4, e5))))

            expect_vector(suppressWarnings(
              errorMetrics(yna, pna, model= "binary", type = "brier")))

            rm(e1, e2, e3, e4, e5, y, p, mm, f1)

})



