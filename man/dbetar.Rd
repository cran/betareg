\name{betar}
\alias{dbetar}
\alias{pbetar}
\alias{qbetar}
\alias{rbetar}

\title{The Beta Distribution in Regression Parameterization}

\description{
  Density, distribution function, quantile function, and random generation
  for the beta distribution in regression parameterization.
}
\usage{
dbetar(x, mu, phi, log = FALSE)

pbetar(q, mu, phi, lower.tail = TRUE, log.p = FALSE)

qbetar(p, mu, phi, lower.tail = TRUE, log.p = FALSE)

rbetar(n, mu, phi)
}
\arguments{
  \item{x, q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mu}{numeric. The mean of the beta distribution.}
  \item{phi}{numeric. The precision parameter of the beta distribution.}
  \item{log, log.p}{logical. If TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If TRUE (default), probabilities are P[X <= x]
    otherwise, P[X > x].}
}

\details{
This is the reparameterization of the beta distribution with mean \code{mu}
and precision \code{phi}, as employed in beta regression. The classic
parameterization of the beta distribution is obtained by setting
\code{shape1 = mu * phi} and \code{shape2 = (1 - mu) * phi},
respectively.
}

\value{
 \code{dbetar} gives the density, \code{pbetar} gives the distribution
     function, \code{qbetar} gives the quantile function, and \code{rbetar}
     generates random deviates.
}

\seealso{\code{\link{dbeta}}, \code{\link{BetaR}}}

\keyword{distribution}
