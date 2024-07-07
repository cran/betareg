\name{BetaR}
\alias{BetaR}
\alias{mean.BetaR}
\alias{variance.BetaR}
\alias{skewness.BetaR}
\alias{kurtosis.BetaR}
\alias{pdf.BetaR}
\alias{log_pdf.BetaR}
\alias{cdf.BetaR}
\alias{quantile.BetaR}
\alias{random.BetaR}
\alias{support.BetaR}
\alias{is_discrete.BetaR}
\alias{is_continuous.BetaR}

\title{Create a Beta Regression Distribution}

\description{
Class and methods for beta distributions in regression specification
using the workflow from the \pkg{distributions3} package.
}

\usage{
BetaR(mu, phi)
}

\arguments{
\item{mu}{numeric. The mean of the beta distribution.}
\item{phi}{numeric. The precision parameter of the beta distribution.}
}

\details{
Alternative parameterization of the classic beta distribution in
terms of its mean \code{mu} and precision parameter \code{phi}.
Thus, the distribution provided by \code{BetaR} is equivalent to
the \code{\link[distributions3]{Beta}} distribution with parameters
\code{alpha = mu * phi} and \code{beta = (1 - mu) * phi}.
}

\value{
A \code{BetaR} distribution object.
}

\seealso{\code{\link{dbetar}}, \code{\link[distributions3]{Beta}}}

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
X <- BetaR(
  mu  = c(0.25, 0.50, 0.75),
  phi = c(1, 1, 2)
)
\donttest{
X

## compute moments of the distribution
mean(X)
variance(X)
skewness(X)
kurtosis(X)

## support interval (minimum and maximum)
support(X)

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
