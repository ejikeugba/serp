library(testthat)
library(serp)

context("check model weights and the errorMetrics function")
wine <- serp::wine

## checks on weight argument
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

})

## check if trace works
expect_output(serp(rating ~ temp + contact, link = "logit",
                   slope = "unparallel", reverse=FALSE,
                   control= list(trace=1),
                   data=wine))

expect_output(serp(rating ~ temp + contact, link = "logit",
                   slope = "unparallel", reverse=FALSE,
                   control= list(trace=2),
                   data=wine))

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

mm <- glm(sample(c(0,1), 50, replace = TRUE) ~ rnorm(50))
expect_error(
  errorMetrics(mm$y, mm$fitted.values, model= "binary", type = "brier"),
  "'actual' must be a factor")

y <- factor(mm$y)
p <- mm$fitted.values
q <- c(NA, p[-1L])

e4 <- suppressWarnings(errorMetrics(y, q, model= "binary", type = "brier"))
e5 <- errorMetrics(y, p, model= "binary", type = "logloss")
e6 <- errorMetrics(y, p, model= "binary", type = "misclass")

expect_false(any(is.na(c(e1, e2, e3, e4, e5, e6))))

rm(e1, e2, e3, e4, e5, e6, y, p, q, mm, f1)

