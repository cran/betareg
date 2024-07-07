\name{XBeta}
\alias{XBeta}
\alias{mean.XBeta}
\alias{variance.XBeta}
\alias{skewness.XBeta}
\alias{kurtosis.XBeta}
\alias{pdf.XBeta}
\alias{log_pdf.XBeta}
\alias{cdf.XBeta}
\alias{quantile.XBeta}
\alias{random.XBeta}
\alias{support.XBeta}
\alias{is_discrete.XBeta}
\alias{is_continuous.XBeta}

\title{Create an Extended-Support Beta Distribution}

\description{
Class and methods for extended-support beta distributions
using the workflow from the \pkg{distributions3} package.
}

\usage{
XBeta(mu, phi, nu = 0)
}

\arguments{
  \item{mu}{numeric. The mean of the underlying beta distribution on [-nu, 1 + nu].}
  \item{phi}{numeric. The precision parameter of the underlying beta
    distribution on [-nu, 1 + nu].}
  \item{nu}{numeric. Exceedence parameter for the support of the underlying
    beta distribution on [-nu, 1 + nu] that is censored to [0, 1].}
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
A \code{XBeta} distribution object.
}

\seealso{\code{\link{dxbeta}}, \code{\link{BetaR}}}

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
X <- XBeta(
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
