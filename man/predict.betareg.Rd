\name{predict.betareg}
\alias{predict.betareg}

\title{Prediction Method for betareg Objects}

\description{
  Extract various types of predictions from beta regression models:
  either on the scale of responses in (0, 1) or the scale of
  the linear predictor.
}

\usage{
\method{predict}{betareg}(object, newdata = NULL,
  type = c("response", "link"), na.action = na.pass, \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict. If omitted, the original observations are used.}
  \item{type}{character indicating type of residuals.}
  \item{na.action}{function determining what should be done with missing values
    in \code{newdata}. The default is to predict \code{NA}.}
  \item{\dots}{currently not used.}
}

\examples{
data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

cbind(
  predict(gy, type = "response"),
  predict(gy, type = "link")
)
}

\keyword{regression}
