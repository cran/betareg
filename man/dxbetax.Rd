\name{xbetax}
\alias{dxbetax}
\alias{pxbetax}
\alias{qxbetax}
\alias{rxbetax}

\title{The Extended-Support Beta Mixture Distribution}

\description{
  Density, distribution function, quantile function, and random generation
  for the extended-support beta mixture distribution (in regression parameterization)
  on [0, 1].
}
\usage{
dxbetax(x, mu, phi, nu = 0, log = FALSE, quad = 20)

pxbetax(q, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20)

qxbetax(p, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20,
  tol = .Machine$double.eps^0.7)

rxbetax(n, mu, phi, nu = 0)
}
\arguments{
  \item{x, q}{numeric. Vector of quantiles.}
  \item{p}{numeric. Vector of probabilities.}
  \item{n}{numeric. Number of observations. If \code{length(n) > 1}, the length is
  taken to be the number required.}
  \item{mu}{numeric. The mean of the underlying beta distribution on [-nu, 1 + nu].}
  \item{phi}{numeric. The precision parameter of the underlying beta
    distribution on [-nu, 1 + nu].}
  \item{nu}{numeric. Mean of the exponentially-distributed exceedence parameter
    for the underlying beta distribution on [-nu, 1 + nu] that is censored to [0, 1].}
  \item{log, log.p}{logical. If TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If TRUE (default), probabilities are P[X <= x]
    otherwise, P[X > x].}
  \item{quad}{numeric. The number of quadrature points for numeric
    integration of the continuous mixture. Alternatively, a matrix with nodes
    and weights for the quadrature points can be specified.}
  \item{tol}{numeric. Accuracy (convergence tolerance) for numerically
    determining quantiles based on \code{\link[stats]{uniroot}} and \code{pxbetax}.}
}

\details{
The extended-support beta mixture distribution is a continuous mixture of
extended-support beta distributions on [0, 1] where the underlying exceedence
parameter is exponentially distributed with mean \code{nu}. Thus, if \code{nu > 0},
the resulting distribution has point masses on the boundaries 0 and 1 with larger
values of \code{nu} leading to higher boundary probabilities. For \code{nu = 0}
(the default), the distribution reduces to the classic beta distribution (in
regression parameterization) without boundary observations.
}

\value{
 \code{dxbetax} gives the density, \code{pxbetax} gives the distribution
     function, \code{qxbetax} gives the quantile function, and \code{rxbetax}
     generates random deviates.
}

\seealso{\code{\link{dxbeta}}, \code{\link{XBetaX}}}

\keyword{distribution}
