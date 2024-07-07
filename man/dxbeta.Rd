\name{xbeta}
\alias{dxbeta}
\alias{pxbeta}
\alias{qxbeta}
\alias{rxbeta}

\title{The Extended-Support Beta Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the extended-support beta distribution (in regression parameterization)
  on [0, 1].
}
\usage{
dxbeta(x, mu, phi, nu = 0, log = FALSE)

pxbeta(q, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE)

qxbeta(p, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE)

rxbeta(n, mu, phi, nu = 0)
}
\arguments{
  \item{x, q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mu}{numeric. The mean of the underlying beta distribution on [-nu, 1 + nu].}
  \item{phi}{numeric. The precision parameter of the underlying beta
    distribution on [-nu, 1 + nu].}
  \item{nu}{numeric. Exceedence parameter for the support of the underlying
    beta distribution on [-nu, 1 + nu] that is censored to [0, 1].}
  \item{log, log.p}{logical. If TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If TRUE (default), probabilities are P[X <= x]
    otherwise, P[X > x].}
}

\details{
In order to obtain an extended-support beta distribution on [0, 1]
an additional exceedence parameter \code{nu} is introduced. If \code{nu > 0},
this scales the underlying beta distribution to the interval [-nu, 1 + nu]
where the tails are subsequently censored to the unit interval [0, 1] with
point masses on the boundaries 0 and 1. Thus, \code{nu} controls how likely
boundary observations are and for \code{nu = 0} (the default), the distribution
reduces to the classic beta distribution (in regression parameterization)
without boundary observations.
}

\value{
 \code{dxbeta} gives the density, \code{pxbeta} gives the distribution
     function, \code{qxbeta} gives the quantile function, and \code{rxbeta}
     generates random deviates.
}

\seealso{\code{\link{dbetar}}, \code{\link{XBeta}}}

\keyword{distribution}
