\name{predict.betareg}
\alias{predict.betareg}
\alias{pit.betareg}
\alias{rootogram.betareg}

\title{Prediction Method for betareg Objects}

\description{
  Extract various types of predictions from beta regression models:
  First, GLM-style predictions on the scale of responses in (0, 1)
  or the scale of the linear predictor are provided. Second, various
  quantities based on the predicted beta distributions are available,
  e.g., moments, quantiles, probabilities, densities, etc.
}

\usage{
\method{predict}{betareg}(object, newdata = NULL,
  type = c("response", "link", "precision", "variance", "parameters",
    "distribution", "density", "probability", "quantile"),
  na.action = na.pass, at = 0.5, elementwise = NULL, \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character indicating type of predictions: fitted means of the response
    (default, \code{"response"} or equivalently \code{"mean"}),
    corresponding linear predictor (\code{"link"}), fitted precision parameter
    phi (\code{"precision"}), fitted variances of the response (\code{"variance"}),
    all \code{"parameters"} of the response distribution, or the corresponding
    \code{"distribution"} object (using the infrastructure from \pkg{distributions3}).
    Finally, standard functions for the distribution can be evaluated (at
    argument \code{at}, see below), namely
    the \code{"density"} (or equivalently \code{"pdf"}),
    the \code{"quantile"} function,
    or the cumulative \code{"probability"} (or equivalently \code{"cdf"}).}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{at}{numeric vector at which the predictions should be evaluated if
    \code{type} specifies a function that takes an additional argument.}
  \item{elementwise}{logical. Should each element of the distribution only be evaluated at the
    corresponding element of \code{at} (\code{elementwise = TRUE}) or at all elements
    in \code{at} (\code{elementwise = FALSE}). Elementwise evaluation is only possible
    if the number of observations is the same as the length of \code{at} and in that case a vector of
    the same length is returned. Otherwise a matrix is returned. The default is to use
    \code{elementwise = TRUE} if possible, and otherwise \code{elementwise = FALSE}.}
  \item{\dots}{further arguments when \code{type} specifies a function, e.g.,
    \code{type = "density"} with \code{log = TRUE} computes log-densities
    (or log-likelihoods).}
}

\details{
  Each prediction for a \code{\link{betareg}} model internally first computes the parameters
  for all observations in \code{newdata} (or the original data if \code{newdata} is missing).
  These parameters correspond to a predicted beta distribution (or extended-support beta distribution)
  for each observation. Then the actual predictions can also be moments of the distributions
  or standard quantities such as densities, cumulative probabilities, or quantiles.
  The latter are computed with the \code{d}/\code{p}/\code{q} functions such as
  \code{\link{dbetar}} (or \code{\link{dxbetax}} or \code{\link{dxbeta}}).
}

\value{
  Either a vector or matrix of predictions with the same number of observations as rows
  in \code{newdata}.
}

\examples{
options(digits = 4)

data("GasolineYield", package = "betareg")

gy2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)

cbind(
  predict(gy2, type = "response"),
  predict(gy2, type = "link"),
  predict(gy2, type = "precision"),
  predict(gy2, type = "variance"),
  predict(gy2, type = "quantile", at = c(0.25, 0.5, 0.75))
)

## evaluate cumulative _p_robabilities for (small) new data set
gyd <- GasolineYield[c(1, 5, 10), ]

## CDF at 0.1 for each observation
predict(gy2, newdata = gyd, type = "probability", at = 0.1)

## CDF at each combination of 0.1/0.2 and observations
predict(gy2, newdata = gyd, type = "probability", at = c(0.1, 0.2))

## CDF at elementwise combinations of 0.1/0.2/0.3 and observations
predict(gy2, newdata = gyd, type = "probability", at = c(0.1, 0.2, 0.3))
predict(gy2, newdata = gyd, type = "probability", at = c(0.1, 0.2, 0.3), elementwise = TRUE)

## CDF at all combinations of 0.1/0.2/0.3 and observations
predict(gy2, newdata = gyd, type = "probability", at = c(0.1, 0.2, 0.3), elementwise = FALSE)
}

\keyword{regression}
