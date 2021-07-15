is_VGAM_installed <- require(VGAM)
library(testthat)
library(serp)

context("unpenalized - compare results from serp with the VGAM::vglm function")
wine <- serp::wine

test_that("cumulative model via serp matches with vglm",
          {
            if (!is_VGAM_installed) skip("VGAM not installed")
            tol <- 1e-06

            #1# parallel slope with logit link
            sp <- serp(rating ~ temp + contact, slope = "parallel",
                       link = "logit", reverse=FALSE, data=wine)
            vm <- vglm(rating ~ temp + contact,
                       family= cumulative(link="logitlink", parallel=TRUE,
                                          reverse=FALSE),data = wine)
            pred.sp <- as.matrix(predict(sp, type="response"))
            pred.vm <- as.matrix(predict(vm, type="response"))

            expect_vector(predict(sp, type="class"))
            expect_equal(coef(sp), coef(vm), check.attributes=FALSE,
                         tolerance=tol)
            expect_equal(pred.sp, pred.vm, check.attributes=FALSE,
                         tolerance=tol)

            expect_error(anova.serp(sp, sp),
                         "duplicate object names are not allowed")

            rm(sp, vm, pred.sp, pred.vm)

            #2# unparallel slope with probit link
            sp <- serp(rating ~ temp + contact, slope = "unparallel",
                       link = "probit", reverse=FALSE, data = wine)
            vm <- vglm(rating ~ temp + contact,
                       family= cumulative(link="probitlink", parallel=FALSE,
                                          reverse=FALSE), data = wine)
            fitted.sp <- as.matrix(sp$fitted.values)
            fitted.vm <- as.matrix(vm@fitted.values)
            expect_equal(fitted.sp, fitted.vm, check.attributes=FALSE,
                         tolerance=tol)
            rm(sp, vm, fitted.sp, fitted.vm)

            #3# parallel slope with reversed cloglog link
            sp <- serp(rating ~ temp + contact, slope = "parallel",
                       link = "cloglog", reverse=TRUE,
                       data = wine)
            vm <- vglm(rating ~ temp + contact,
                       family=cumulative(link="clogloglink", parallel=TRUE,
                                         reverse=TRUE), data = wine)
            loglik.sp <- sp$logLik
            loglik.vm <- logLik(vm)
            expect_equal(loglik.sp, loglik.vm, check.attributes=FALSE,
                         tolerance=tol)
            rm(sp, vm, loglik.sp, loglik.vm)

            #4# unparallel slope with reverse logit link
            sp <- serp(rating ~ temp + contact, slope = "unparallel",
                       link = "logit", reverse=TRUE, data = wine)
            vm <- vglm(rating ~ temp + contact,
                       family=cumulative(link="logitlink", parallel=FALSE,
                                         reverse=TRUE), data = wine)
            iter.sp <- sp$iter
            iter.vm <- vm@iter
            expect_equal(iter.sp, iter.vm, check.attributes=FALSE,
                         tolerance=tol)
            rm(sp, vm, iter.sp, iter.vm)

            #5# partial or semi-parallel slope with cauchit link
            sp <- serp(rating ~ temp * contact, slope = "partial",
                       link = "cauchit", globalEff= ~temp, data = wine)
            vm <- vglm(rating ~ temp * contact,
                       family=cumulative(link="cauchitlink", parallel=FALSE ~contact),
                       data = wine)
            dev.sp <- deviance(sp)
            dev.vm <- deviance(vm)

            cof.sp <- summary(sp)$coefficient[,1L]
            cof.vm <- coef(vm)
            expect_equal(cof.sp, cof.vm, check.attributes=FALSE,
                         tolerance=1e-06)

            expect_error(
              serp(rating ~ temp * contact, slope = "partial",
                   link = "cauchit", globalEff= "temp", data = wine),
              "no object of class formula used in globalEff")

            expect_error(
              serp(rating ~ temp * contact, slope = "partial",
                   link = "cauchit", data = wine),
              "'globalEff' is unspecified")


            rm(sp, vm, dev.sp, dev.vm, cof.sp, cof.vm, tol)
          })





