\name{Beta4}
\alias{Beta4}
\alias{mean.Beta4}
\alias{variance.Beta4}
\alias{skewness.Beta4}
\alias{kurtosis.Beta4}
\alias{pdf.Beta4}
\alias{log_pdf.Beta4}
\alias{cdf.Beta4}
\alias{quantile.Beta4}
\alias{random.Beta4}
\alias{support.Beta4}
\alias{is_discrete.Beta4}
\alias{is_continuous.Beta4}

\title{Create a 4-Parameter Beta Distribution}

\description{
Class and methods for 4-parameter beta distributions in regression specification
using the workflow from the \pkg{distributions3} package.
}

\usage{
Beta4(mu, phi, theta1 = 0, theta2 = 1 - theta1)
}

\arguments{
  \item{mu}{numeric. The mean of the beta distribution that is extended to
    support [theta1, theta2].}
  \item{phi}{numeric. The precision parameter of the beta distribution that is
    extended to support [theta1, theta2].}
  \item{theta1, theta2}{numeric. The minimum and maximum, respectively,
    of the 4-parameter beta distribution. By default a symmetric support is
    chosen by \code{theta2 = 1 - theta1} which reduces to the classic
    beta distribution because of the default \code{theta1 = 0}.}
}

\details{
The distribution is obtained by a linear transformation of a beta-distributed
random variable with intercept \code{theta1} and slope \code{theta2 - theta1}.
}

\value{
A \code{Beta4} distribution object.
}

\seealso{\code{\link{dbeta4}}, \code{\link{BetaR}}}

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
X <- Beta4(
  mu  = c(0.25, 0.50, 0.75),
  phi = c(1, 1, 2),
  theta1 = c(0, -0.1, -0.1),
  theta2 = c(1, 1.1, 1.5)
)
\donttest{
X

## compute moments of the distribution
mean(X)
variance(X)

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

## cdf() and quantile() are inverses
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
