\name{gleverage}
\alias{gleverage}
\alias{gleverage.betareg}

\title{Generalized Leverage Values}

\description{
  Compute the generalized leverages values for fitted models.
}

\usage{
gleverage(model, \dots)
}

\arguments{
  \item{model}{a model object.}
  \item{\dots}{further arguments passed to methods.}
}

\value{
\code{gleverage} is a new generic for computing generalized leverage values as suggested by
Wei, Hu, and Fung (1998). Currently, there is only a method for \code{betareg} models, implementing
the formulas from Rocha and Simas (2010) which are consistent with the formulas from
Ferrari and Cribari-Neto (2004) for the fixed dispersion case.

Currently, the vector of generalized leverages requires computations and
storage of order \eqn{n \times n}{n x n}.
}

\references{
Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Rocha, A.V., and Simas,  A.B. (2010).
  Influence Diagnostics in a General Class of Beta Regression Models.
  \emph{Test}, Forthcoming. \url{http://dx.doi.org/10.1007/s11749-010-0189-z}

Wei, B.-C., and Hu, Y.-Q., and Fung, W.-K. (1998).
  Generalized Leverage and Its Applications.
  \emph{Scandinavian Journal of Statistics}, \bold{25}, 25--37.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")
gy <- betareg(yield ~ batch + temp, data = GasolineYield)
gleverage(gy)
}

\keyword{regression}
