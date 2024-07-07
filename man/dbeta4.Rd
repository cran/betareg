\name{beta4}
\alias{dbeta4}
\alias{pbeta4}
\alias{qbeta4}
\alias{rbeta4}

\title{The 4-Parameter Beta Distribution in Regression Parameterization}

\description{
  Density, distribution function, quantile function, and random generation
  for the 4-parameter beta distribution in regression parameterization.
}
\usage{
dbeta4(x, mu, phi, theta1 = 0, theta2 = 1 - theta1, log = FALSE)

pbeta4(q, mu, phi, theta1 = 0, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE)

qbeta4(p, mu, phi, theta1 = 0, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE)

rbeta4(n, mu, phi, theta1 = 0, theta2 = 1 - theta1)
}
\arguments{
  \item{x, q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mu}{numeric. The mean of the beta distribution that is extended to
    support [theta1, theta2].}
  \item{phi}{numeric. The precision parameter of the beta distribution that is
    extended to support [theta1, theta2].}
  \item{theta1, theta2}{numeric. The minimum and maximum, respectively,
    of the 4-parameter beta distribution. By default a symmetric support is
    chosen by \code{theta2 = 1 - theta1} which reduces to the classic
    beta distribution because of the default \code{theta1 = 0}.}
  \item{log, log.p}{logical. If TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If TRUE (default), probabilities are P[X <= x]
    otherwise, P[X > x].}
}

\details{
The distribution is obtained by a linear transformation of a beta-distributed
random variable with intercept \code{theta1} and slope \code{theta2 - theta1}.
}

\value{
 \code{dbeta4} gives the density, \code{pbeta4} gives the distribution
     function, \code{qbeta4} gives the quantile function, and \code{rbeta4}
     generates random deviates.
}

\seealso{\code{\link{dbetar}}, \code{\link{Beta4}}}

\keyword{distribution}
