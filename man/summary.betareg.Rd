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
\method{summary}{betareg}(object, phi = NULL, type = "sweighted2", \dots)

\method{coef}{betareg}(object, model = c("full", "mean", "precision"), phi = NULL, \dots)
\method{vcov}{betareg}(object, model = c("full", "mean", "precision"), phi = NULL, \dots)
\method{bread}{betareg}(x, phi = NULL, \dots)
\method{estfun}{betareg}(x, phi = NULL, \dots)
}

\arguments{
  \item{object, x}{fitted model object of class \code{"betareg"}.}
  \item{phi}{logical indicating whether the parameters in the precision model
    (for phi) should be reported as full model parameters (\code{TRUE}) or
    nuisance parameters (\code{FALSE}). The default is taken from
    \code{object$phi}.}
  \item{type}{character specifying type of residuals to be included in the
    summary output, see \code{\link{residuals.betareg}}.}
  \item{model}{character specifying for which component of the model coefficients/covariance
    should be extracted. (Only used if \code{phi} is \code{NULL}.)}
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
  using the associated \code{print} method. Note that the default residuals
  \code{"sweighted2"} might be burdensome to compute in large samples and hence might
  need modification in such applications.

  A \code{\link[stats]{logLik}} method is provided, hence \code{\link[stats]{AIC}}
  can be called to compute information criteria.
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Simas, A.B., and Barreto-Souza, W., and Rocha, A.V. (2010).
  Improved Estimators for a General Class of Beta Regression Models.
  \emph{Computational Statistics & Data Analysis}, \bold{54}(2), 348--366.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)

summary(gy2)
coef(gy2)
vcov(gy2)
logLik(gy2)
AIC(gy2)

coef(gy2, model = "mean")
coef(gy2, model = "precision")
summary(gy2, phi = FALSE)
}

\keyword{regression}
