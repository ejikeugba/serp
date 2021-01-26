library(testthat)
library(serp)


context("Penalized - checks if serp works properly on cumulative models")

## To avoid object
## masking
########################
wine <- serp::wine
########################
test_that("estreem shrinkage with serp results to the
           proportional odds model (POM) while zero
           shrinkage returns the non-proportional odds
           model (NPOM)",
{
  tol <- 1e-06
  #1# From NPOM to POM
  sp1 <- serp(rating ~ temp + contact, link = "logit",
              slope = "penalize", tuneMethod = "user",
              lambda = 1e10, reverse=FALSE, data=wine)
  cof1 <- unique(round(as.numeric(coef(sp1)), 9L))

  sp2 <- serp(rating ~ temp + contact, link = "logit",
              slope = "parallel", reverse=FALSE, data = wine)
  cof2 <- coef(sp2)
  expect_equal(cof1, cof2, check.attributes=FALSE,
               tolerance=tol)
  rm(sp1, sp2, cof1, cof2)

  #2# Zero shrinkage yields NPOM
  sp1 <- serp(rating ~ temp + contact, link = "logit",
              slope = "penalize", tuneMethod = "user",
              lambda = 0, reverse=FALSE, data=wine)
  sp2 <- serp(rating ~ temp + contact, link = "logit",
              slope = "unparallel", reverse=FALSE,
              data=wine)

  #3# checks on anova, confinct and vcov functions

  expect_equal(coef(sp1), coef(sp2), check.attributes=FALSE,
               tolerance=tol)

  expect_error(serp:::anova.serp(sp1),
               "no anova implementation yet for a single 'serp' object")
  expect_false(inherits(try(serp:::anova.serp(sp1, sp2)), 'try-error'))


  expect_error(serp:::confint.serp(sp1, sp2), "one object at a time allowed")

  expect_false(inherits(try(serp:::confint.serp(sp1)), 'try-error'))

  expect_error(serp:::vcov.serp(sp1, sp2), "one object at a time allowed")

  expect_false(inherits(try(serp:::vcov.serp(sp1)), 'try-error'))


  expect_output(serp:::print.serp(sp1))
  expect_output(serp:::print.summary.serp(serp:::summary.serp(sp1)))

  rm(sp1, sp2)
})


##
test_that("serp returns expected results for different tuneMethod",
      {
 expect_vector(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       lambdaGrid = c(1,2),
       tuneMethod = "deviance",
       gridType = "discrete",
       data = wine)$lambda
)

 expect_vector(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       lambdaGrid = c(1,2),
       tuneMethod = "cv",
       data = wine,
  )$slope)

 expect_vector(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       tuneMethod = "finite",
       data = wine,
  )$value)

 expect_vector(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = TRUE,
       globalEff = ~ temp,
       tuneMethod = "finite",
       data = wine,
  )$value)

 expect_error(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       tuneMethod = "deviance",
       lambdaGrid = c(0,1), control = list(trace=-1),
       data = wine),
  "maxits, eps, maxpen, minP, maxAdjIter, max.half.iter and relTol should all be numeric and non-negative")
})



## checks on lambda and lambdaGrid
test_that("lambda is a single numeric and non-negative value and
           that lambdaGrid contains the right input values",
    {
 expect_error(
  serp(rating ~ temp + contact, slope = "penalize",
       link = "loglog", reverse = TRUE, tuneMethod = "user",
       lambda = c(0.3,5), data = wine),
  "lambda should be a single numeric and non-negative value")

 expect_error(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       tuneMethod = "deviance",
       lambdaGrid = c(-3,4),
       data = wine),
  "lambdaGrid must be a non-negative numeric vector of length > 1")
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


## checks on predict function
pred1 <- predict(
  serp(rating ~ temp + contact, link = "logit",
       slope = "penalize", reverse = FALSE,
       globalEff = ~ temp + contact, data = wine),
       type = 'response', newdata=head(wine))

pred2 <- predict(
  serp(rating ~ temp + contact, link = "logit",
       slope = "parallel", reverse = FALSE,
       data = wine), type = 'response', newdata=head(wine))
expect_equal(pred1, pred2)

