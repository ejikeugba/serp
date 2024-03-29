#' Bitterness of wine dataset
#'
#' The \code{wine} dataset adopted from Randall(1989),
#' represents the outcome of a factorial experiment on factors
#' determining the bitterness of wine. Two treatment factors
#' (temperature and contact) with two levels each are provided,
#' with the rating of wine taken on a continuous scale in the interval
#' from 0 (none) to 100 (intense). These were subsequently grouped
#' into five ordered categories ranging from 1 = 'least bitter'
#' to 5 = 'most bitter'. Altogether, nine different judges assessed
#' wine from two bottles and out of the four treatment conditions,
#' making a total of 72 observations.
#' @docType data
#' @usage wine
#' @format A data frame with 72 rows and 6 variables:
#' @return \item{\code{response}}{
#'    scorings of wine bitterness on a 0---100 continuous scale.
#'  }
#' @return \item{\code{rating}}{
#'    ordered factor with 5 levels; a grouped version of \code{response}.
#'  }
#' @return \item{\code{contact}}{
#'    factor with two levels (\code{"no"} and \code{"yes"}).
#'  }
#' @return \item{\code{temp}}{
#'    temperature: factor with two levels.
#'  }
#' @return \item{\code{judge}}{
#'    factor with nine levels.
#'  }
#' @return \item{\code{bottle}}{
#'    factor with eight levels.
#'  }
#' @keywords dataset
#' @source Taken from Randall (1989).
#' @references
#' Randall, J (1989). The analysis of sensory data by generalized linear
#'     model. \emph{Biometrical Journal}, 31, 781--793.
#'     https://doi.org/10.1002/bimj.4710310703
#'
#' @examples
#'
#' \dontrun{
#' str(wine)
#' head(wine)
#' }
#'
"wine"
