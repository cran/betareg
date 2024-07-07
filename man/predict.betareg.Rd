\name{predict.betareg}
\alias{predict.betareg}
\alias{pit.betareg}
\alias{rootogram.betareg}

\title{Prediction Method for betareg Objects}

\description{
  Extract various types of predictions from beta regression models:
  either on the scale of responses in (0, 1) or the scale of
  the linear predictor.
}

\usage{
\method{predict}{betareg}(object, newdata = NULL,
  type = c("response", "link", "precision", "variance", "parameters",
    "density", "probability", "quantile"),
  na.action = na.pass, at = 0.5, \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character indicating type of predictions: fitted means of response (\code{"response"}),
    corresponding linear predictor (\code{"link"}), fitted precision parameter
    phi (\code{"precision"}), fitted variances of response (\code{"variance"}),
    or fitted quantile(s) of the response distribution (\code{"quantile"}).}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{at}{numeric vector indicating the level(s) at which quantiles
    should be predicted (only if \code{type = "quantile"}), defaulting
    to the median \code{at = 0.5}.}
  \item{\dots}{currently not used.}
}

\details{
  FIXME: Update to extended type and at processing.
  
  FIXME: Add comments about pit and rootogram.
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
## CDF at pairwise combinations of 0.1/0.2/0.3 and observations
predict(gy2, newdata = gyd, type = "probability", at = c(0.1, 0.2, 0.3))
## CDF at all combinations of 0.1/0.2/0.3 and observations
predict(gy2, newdata = gyd, type = "probability", at = rbind(c(0.1, 0.2, 0.3)))
}

\keyword{regression}
