\name{Beta01}
\alias{Beta01}
\alias{mean.Beta01}
\alias{variance.Beta01}
\alias{skewness.Beta01}
\alias{kurtosis.Beta01}
\alias{pdf.Beta01}
\alias{log_pdf.Beta01}
\alias{cdf.Beta01}
\alias{quantile.Beta01}
\alias{random.Beta01}
\alias{support.Beta01}
\alias{is_discrete.Beta01}
\alias{is_continuous.Beta01}

\title{Create a Zero- and/or One-Inflated Beta Distribution}

\description{
Class and methods for zero- and/or one-inflated beta distributions in regression specification
using the workflow from the \pkg{distributions3} package.
}

\usage{
Beta01(mu, phi, p0 = 0, p1 = 0)
}

\arguments{
  \item{mu}{numeric. The mean of the beta distribution (on the open unit interval).}
  \item{phi}{numeric. The precision parameter of the beta distribution.}
  \item{p0}{numeric. The probability for an observation of zero (often referred
    to as zero inflation).}
  \item{p1}{numeric. The probability for an observation of one (often referred
    to as one inflation).}
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
A \code{Beta01} distribution object.
}

\seealso{\code{\link{dbeta01}}, \code{\link{BetaR}}}

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
X <- Beta01(
  mu  = c(0.25, 0.50, 0.75),
  phi = c(1, 1, 2),
  p0 = c(0.1, 0, 0),
  p1 = c(0, 0, 0.3)
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

## point mass probabilities (if any) on boundary
cdf(X, 0, lower.tail = TRUE)
cdf(X, 1, lower.tail = FALSE)

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
