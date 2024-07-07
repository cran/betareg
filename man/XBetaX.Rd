\name{XBetaX}
\alias{XBetaX}
\alias{mean.XBetaX}
\alias{variance.XBetaX}
\alias{skewness.XBetaX}
\alias{kurtosis.XBetaX}
\alias{pdf.XBetaX}
\alias{log_pdf.XBetaX}
\alias{cdf.XBetaX}
\alias{quantile.XBetaX}
\alias{random.XBetaX}
\alias{support.XBetaX}
\alias{is_discrete.XBetaX}
\alias{is_continuous.XBetaX}

\title{Create an Extended-Support Beta Mixture Distribution}

\description{
Class and methods for extended-support beta distributions
using the workflow from the \pkg{distributions3} package.
}

\usage{
XBetaX(mu, phi, nu = 0)
}

\arguments{
  \item{mu}{numeric. The mean of the underlying beta distribution on [-nu, 1 + nu].}
  \item{phi}{numeric. The precision parameter of the underlying beta
    distribution on [-nu, 1 + nu].}
  \item{nu}{numeric. Mean of the exponentially-distributed exceedence parameter
    for the underlying beta distribution on [-nu, 1 + nu] that is censored to [0, 1].}
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
A \code{XBetaX} distribution object.
}

\seealso{\code{\link{dxbetax}}, \code{\link{XBeta}}}

\examples{
\dontshow{ if(!requireNamespace("distributions3")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and random seed
library("distributions3")
set.seed(6020)

## three beta distributions
X <- XBetaX(
  mu  = c(0.25, 0.50, 0.75),
  phi = c(1, 1, 2),
  nu = c(0, 0.1, 0.2)
)
\donttest{
X

## compute moments of the distribution
mean(X)
variance(X)

## support interval (minimum and maximum)
support(X)

## it is only continuous when there are no point masses on the boundary
is_continuous(X)
cdf(X, 0)
cdf(X, 1, lower.tail = FALSE)

## simulate random variables
random(X, 5)

## histograms of 1,000 simulated observations
x <- random(X, 1000)
hist(x[1, ])
hist(x[2, ])
hist(x[3, ])

## probability density function (PDF) and log-density (or log-likelihood)
x <- c(0.25, 0.5, 0.75)
pdf(X, x)
pdf(X, x, log = TRUE)
log_pdf(X, x)

## cumulative distribution function (CDF)
cdf(X, x)

## quantiles
quantile(X, 0.5)

## cdf() and quantile() are inverses (except at censoring points)
cdf(X, quantile(X, 0.5))
quantile(X, cdf(X, 1))

## all methods above can either be applied elementwise or for
## all combinations of X and x, if length(X) = length(x),
## also the result can be assured to be a matrix via drop = FALSE
p <- c(0.05, 0.5, 0.95)
quantile(X, p, elementwise = FALSE)
quantile(X, p, elementwise = TRUE)
quantile(X, p, elementwise = TRUE, drop = FALSE)

## compare theoretical and empirical mean from 1,000 simulated observations
cbind(
  "theoretical" = mean(X),
  "empirical" = rowMeans(random(X, 1000))
)
}
}
