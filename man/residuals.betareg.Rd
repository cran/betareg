\name{residuals.betareg}
\alias{residuals.betareg}

\title{Residuals Method for betareg Objects}

\description{
  Extract various types of residuals from beta regression models:
  raw response residuals (observed - fitted), Pearson residuals (raw residuals scaled by
  square root of variance function), deviance residuals (scaled log-likelihood contributions),
  and different kinds of weighted residuals suggested by Espinheira et al. (2008).
}

\usage{
\method{residuals}{betareg}(object,
  type = c("sweighted2", "deviance", "pearson", "response", "weighted", "sweighted"),
  \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{type}{character indicating type of residuals.}
  \item{\dots}{currently not used.}
}

\details{
  The definitions of all residuals are provided in Espinheira et al. (2008):
  Equation 2 for \code{"pearson"}, last equation on page 409 for \code{"deviance"},
  Equation 6 for \code{"weighted"}, Equation 7 for \code{"sweighted"}, and
  Equation 8 for \code{"sweighted2"}.
  
  Espinheira et al. (2008) recommend to use \code{"sweighted2"}, hence this is 
  the default in the \code{residuals()} method.
}

\references{
Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008).
  On Beta Regression Residuals.
  \emph{Journal of Applied Statistics}, \bold{35}(4), 407--419.

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

gy_res <- cbind(
  residuals(gy, type = "pearson"),
  residuals(gy, type = "deviance"),
  residuals(gy, type = "response"),
  residuals(gy, type = "weighted"),
  residuals(gy, type = "sweighted"),
  residuals(gy, type = "sweighted2")
)
colnames(gy_res) <- c("pearson", "deviance", "response",
  "weighted", "sweighted", "sweighted2")
pairs(gy_res)
}

\keyword{regression}
