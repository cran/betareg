\name{summary.betareg}

\alias{print.betareg}
\alias{summary.betareg}
\alias{print.summary.betareg}

\alias{coef.betareg}
\alias{vcov.betareg}
\alias{bread.betareg}
\alias{estfun.betareg}

\alias{coeftest.betareg}
\alias{logLik.betareg}
\alias{terms.betareg}
\alias{model.frame.betareg}
\alias{model.matrix.betareg}

\alias{cooks.distance.betareg}
\alias{hatvalues.betareg}

\title{Methods for betareg Objects}

\description{
  Methods for extracting information from fitted beta
  regression model objects of class \code{"betareg"}.
}

\usage{
\method{summary}{betareg}(object, phi = NULL, \dots)

\method{coef}{betareg}(object, phi = NULL, \dots)
\method{vcov}{betareg}(object, phi = NULL, \dots)
\method{bread}{betareg}(x, phi = NULL, \dots)
\method{estfun}{betareg}(x, phi = NULL, \dots)
}

\arguments{
  \item{object, x}{fitted model object of class \code{"betareg"}.}
  \item{phi}{logical indicating whether the precision/dispersion parameter
    phi should be reported as a full model parameter (\code{TRUE}) or
    a nuisance parameter (\code{FALSE}). The default is taken from
    \code{object$phi}.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"betareg"}, including methods to the generic functions
  \code{\link[base]{print}} and \code{\link[base]{summary}} which print the estimated
  coefficients along with some further information. The \code{summary} in particular
  supplies partial Wald tests based on the coefficients and the covariance matrix.
  As usual, the \code{summary} method returns an object of class \code{"summary.betareg"}
  containing the relevant summary statistics which can subsequently be printed
  using the associated \code{print} method.

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\references{
Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

summary(gy)
coef(gy)
vcov(gy)
logLik(gy)
AIC(gy)
}

\keyword{regression}
