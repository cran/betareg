\name{CarTask}
\alias{CarTask}
\docType{data}
\title{
  Partition-primed Probability Judgement Task for Car Dealership
}
\description{
  In this study participants were asked to judge how likely it is that a
  customer trades in a coupe or that a customer buys a car form a
  specific salesperson out of four possible salespersons.
}
\usage{data(CarTask)}
\format{
  A data frame with 155 observations on the following 3 variables.
  \describe{
    \item{\code{task}}{a factor with levels \code{Car} and
      \code{Salesperson} indicating the condition.}
    \item{\code{probability}}{a numeric vector of the estimated probability.}
    \item{\code{NFCCscale}}{a numeric vector of the NFCC scale.}
  }
}
\details{
  All participants in the study were undergraduate students at The
  Australian National University, some of whom obtained course credit in
  first-year Psychology for their participation in the study.

  The NFCC scale is a combined scale of the Need for Closure and Need
  for Certainty scales which are strongly correlated.
  
  For \code{task} the questions were:
  \describe{
    \item{Car}{What is the probability that a customer trades in a coupe?}
    \item{Salesperson}{What is the probability that a customer buys a
      car from Carlos?}
  }
}
\source{
  Taken from \url{http://dl.dropbox.com/u/1857674/betareg/betareg.html}.
}
\references{
  Smithson, M., Merkle, E.C., and Verkuilen, J. (in press). Beta
  Regression Finite Mixture Models of Polarization and
  Priming. \emph{Journal of Educational and Behavioral Statistics}.

  Smithson, M., and Segale, C. (2009). Partition Priming in Judgments of
  Imprecise Probabilities. \emph{Journal of Statistical Theory and
  Practice}, \bold{3}(1), 169--181.
}
\examples{
data("CarTask", package = "betareg")
library("flexmix")
car_betamix <- betamix(probability ~ 1, data = CarTask, k = 3,
  extra_components = list(extraComponent(type = "uniform", coef = 1/2,
  delta = 0.01), extraComponent(type = "uniform", coef = 1/4, delta = 0.01)),
  FLXconcomitant = FLXPmultinom(~ task))
}
\keyword{datasets}
