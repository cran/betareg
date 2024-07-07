\name{beta01}
\alias{dbeta01}
\alias{pbeta01}
\alias{qbeta01}
\alias{rbeta01}

\title{The Zero- and/or One-Inflated Beta Distribution in Regression Parameterization}

\description{
  Density, distribution function, quantile function, and random generation
  for the zero- and/or one-inflated beta distribution in regression parameterization.
}
\usage{
dbeta01(x, mu, phi, p0 = 0, p1 = 0, log = FALSE)

pbeta01(q, mu, phi, p0 = 0, p1 = 0, lower.tail = TRUE, log.p = FALSE)

qbeta01(p, mu, phi, p0 = 0, p1 = 0, lower.tail = TRUE, log.p = FALSE)

rbeta01(n, mu, phi, p0 = 0, p1 = 0)
}
\arguments{
  \item{x, q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mu}{numeric. The mean of the beta distribution (on the open unit interval).}
  \item{phi}{numeric. The precision parameter of the beta distribution.}
  \item{p0}{numeric. The probability for an observation of zero (often referred
    to as zero inflation).}
  \item{p1}{numeric. The probability for an observation of one (often referred
    to as one inflation).}
  \item{log, log.p}{logical. If TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If TRUE (default), probabilities are P[X <= x]
    otherwise, P[X > x].}
}

\details{
The zero- and/or one-inflated beta distribution is obtained by adding point
masses at zero and/or one to a standard beta distribution. 

Note that the support of the standard beta distribution is the open unit
interval where values of exactly zero or one cannot occur. Thus, the inflation
jargon is rather misleading as there is no probability that could be inflated.
It is rather a hurdle or two-part (or three-part) model.
}

\value{
 \code{dbeta01} gives the density, \code{pbeta01} gives the distribution
     function, \code{qbeta01} gives the quantile function, and \code{rbeta01}
     generates random deviates.
}

\seealso{\code{\link{dbetar}}, \code{\link{Beta01}}}

\keyword{distribution}
