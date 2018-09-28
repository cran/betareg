\name{ImpreciseTask}
\alias{ImpreciseTask}
\docType{data}
\title{
  Imprecise Probabilities for Sunday Weather and Boeing Stock Task
}
\description{
  In this study participants were asked to estimate upper and lower
  probabilities for event to occur and not to occur.
}
\usage{data(ImpreciseTask)}
\format{
  A data frame with 242 observations on the following 3 variables.
  \describe{
    \item{\code{task}}{a factor with levels \code{Boeing stock} and \code{Sunday weather}.}
    \item{\code{location}}{a numeric vector of the average of the lower
      estimate for the event not to occur and the upper estimate for the
      event to occur.}
    \item{\code{difference}}{a numeric vector of the differences of the
      lower and upper estimate for the event to occur.}
  }
}
\details{
  All participants in the study were either first- or second-year
  undergraduate students in psychology, none of whom had a strong
  background in probability or were familiar with imprecise probability
  theories. 

  For the sunday weather task see \code{\link{WeatherTask}}. For the Boeing
  stock task participants were asked to estimate the probability that
  Boeing's stock would rise more than those in a list of 30 companies.

  For each task participants were asked to provide lower and upper
  estimates for the event to occur and not to occur.  
}

\source{
  Taken from Smithson et al. (2011) supplements.
}

\references{
  Smithson, M., Merkle, E.C., and Verkuilen, J. (2011). Beta
  Regression Finite Mixture Models of Polarization and Priming.
  \emph{Journal of Educational and Behavioral Statistics}, \bold{36}(6), 804--831.
  \doi{10.3102/1076998610396893}

  Smithson, M., and Segale, C. (2009). Partition Priming in Judgments of
  Imprecise Probabilities. \emph{Journal of Statistical Theory and
  Practice}, \bold{3}(1), 169--181.
}
\examples{
data("ImpreciseTask", package = "betareg")
library("flexmix")
wt_betamix <- betamix(location ~ difference * task, data = ImpreciseTask, k = 2,
  extra_components = extraComponent(type = "betareg", coef =
    list(mean = 0, precision = 8)),
  FLXconcomitant = FLXPmultinom(~ task))
}
\keyword{datasets}
