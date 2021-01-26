library(testthat)
library(serp)
context("checks the availability of formula and response in model")
wine <- serp::wine

##
test_that("formula is specified in serp",
{
  expect_error(
    serp(globalEff=~temp, link = "probit", reverse=TRUE,
         slope = "parallel", data = wine),
    "Model needs a formula")
  expect_error(
    serp(data = wine), "Model needs a formula")
  expect_vector(serp(rating~ 1, link = "cauchit", slope = "penalize",
                     tuneMethod = "finite", reverse=TRUE,
                     data = wine)$coef)
})

##
test_that("response exist, is not in predictor or in 'GlobalEff'",
{
  expect_error(
    serp(~ temp + contact, link = "cauchit", slope = "penalize",
         tuneMethod = "finite", reverse=TRUE, data = wine),
    "response missing in formula")
  expect_error(
    serp(rating ~ temp + rating + contact, link = "logit",
         slope = "penalize", reverse=TRUE, data = wine),
    "response not allowed as predictor")
  expect_error(
    serp(rating ~ temp + contact, link = "logit",
         slope = "penalize", reverse=TRUE, tuneMethod = "cv",
         globalEff = ~rating, data = wine),
    "response not allowed in 'globalEff'")
})
