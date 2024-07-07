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
\method{residuals}{betareg}(object, type = c("quantile",
  "deviance", "pearson", "response", "weighted", "sweighted", "sweighted2"),
  \dots)
}

\arguments{
  \item{object}{fitted model object of class \code{"betareg"}.}
  \item{type}{character indicating type of residuals.}
  \item{\dots}{currently not used.}
}

\details{
  The default residuals (starting from version 3.2-0) are quantile residuals
  as proposed by Dunn and Smyth (1996) and explored in the context of beta
  regression by Pereira (2017). In case of extended support beta regression
  with boundary observations at 0 and/or 1, the quantile residuals for the
  boundary observations are randomized.

  The definitions of all other residuals are provided in Espinheira et al. (2008):
  Equation 2 for \code{"pearson"}, last equation on page 409 for \code{"deviance"},
  Equation 6 for \code{"weighted"}, Equation 7 for \code{"sweighted"}, and
  Equation 8 for \code{"sweighted2"}.
  
  Espinheira et al. (2008) recommend to use \code{"sweighted2"}, hence this was 
  the default prior to version 3.2-0. However, these are rather burdensome to
  compute because they require operations of \eqn{O(n^2)} and hence are typically
  prohibitively costly in large sample. Also they are not available for
  extended support beta regression. Finally, Pereira (2017) found quantile
  residuals to have better distributional properties.
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \doi{10.18637/jss.v034.i02}

Dunn, P.K., and Smyth, G.K. (1996). Randomized Quantile Residuals.
  \emph{Journal of Computational and Graphical Statistics}, \bold{5}(3), 236--244.
  \doi{10.2307/1390802}

Espinheira, P.L., Ferrari, S.L.P., and Cribari-Neto, F. (2008).
  On Beta Regression Residuals.
  \emph{Journal of Applied Statistics}, \bold{35}(4), 407--419.
  \doi{10.1080/02664760701834931}

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
  \doi{10.1080/0266476042000214501}

Pereira, G.H.A. (2017). On Quantile Residuals in Beta Regression.
  \emph{Communications in Statistics -- Simulation and Computation}, \bold{48}(1), 302--316.
  \doi{10.1080/03610918.2017.1381740}
}

\seealso{\code{\link{betareg}}}

\examples{
options(digits = 4)

data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

gy_res <- cbind(
  "quantile"   = residuals(gy, type = "quantile"),
  "pearson"    = residuals(gy, type = "pearson"),
  "deviance"   = residuals(gy, type = "deviance"),
  "response"   = residuals(gy, type = "response"),
  "weighted"   = residuals(gy, type = "weighted"),
  "sweighted"  = residuals(gy, type = "sweighted"),
  "sweighted2" = residuals(gy, type = "sweighted2")
)
pairs(gy_res)

cor(gy_res)
}

\keyword{regression}
