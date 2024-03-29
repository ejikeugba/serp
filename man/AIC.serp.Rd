% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serp.method.R
\name{AIC.serp}
\alias{AIC.serp}
\title{AIC for a fitted serp object}
\usage{
\method{AIC}{serp}(object, ..., k = 2)
}
\arguments{
\item{object}{An object of class \code{serp}.}

\item{...}{additional arguments.}

\item{k}{fixed value equal to 2.}
}
\value{
A single numeric value of the model AIC.
}
\description{
Returns the akaike information criterion of a fitted object of class
\code{serp}. For the penalized slope, the effective degrees of freedom (edf)
is obtained from the trace of the generalized hat matrix which depends on
the tuning parameter.
}
\examples{
library(serp)
m <- serp(rating ~ temp + contact, slope = "parallel", link = "probit",
          data = wine)
AIC(m)
}
\seealso{
\code{\link{serp}}, \code{\link{BIC.serp}}, \code{\link{coef.serp}},
\code{\link{logLik.serp}},
}
