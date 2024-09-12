\name{WeatherTask}
\alias{WeatherTask}

\title{Weather Task with Priming and Precise and Imprecise Probabilities}

\description{
  In this study participants were asked to judge how likely Sunday is to
  be the hottest day of the week.
}
\usage{data("WeatherTask", package = "betareg")}
\format{
  A data frame with 345 observations on the following 3 variables.
  \describe{
    \item{\code{priming}}{a factor with levels \code{two-fold} (case
      prime) and \code{seven-fold} (class prime).}
    \item{\code{eliciting}}{a factor with levels \code{precise} and
      \code{imprecise} (lower and upper limit).}
    \item{\code{agreement}}{a numeric vector, probability indicated by
      participants or the average between minimum and maximum
      probability indicated.}
  }
}
\details{
  All participants in the study were either first- or second-year
  undergraduate students in psychology, none of whom had a strong
  background in probability or were familiar with imprecise probability
  theories. 
  
  For \code{priming} the questions were:
  \describe{
    \item{two-fold}{[What is the probability that] the temperature at
      Canberra airport on Sunday will be higher than every other day
      next week?}
    \item{seven-fold}{[What is the probability that] the highest
      temperature of the week at Canberra airport will occur on Sunday?}
  }
  For \code{eliciting} the instructions were if
  \describe{
    \item{precise}{to assign a probability estimate,}
    \item{imprecise}{to assign a lower and upper probability estimate.}
  }
}

\source{
  Taken from Smithson et al. (2011) supplements.
}

\references{
  Smithson M, Merkle EC, Verkuilen J (2011). Beta
  Regression Finite Mixture Models of Polarization and Priming.
  \emph{Journal of Educational and Behavioral Statistics}, \bold{36}(6), 804--831.
  \doi{10.3102/1076998610396893}

  Smithson M, Segale C (2009). Partition Priming in Judgments of
  Imprecise Probabilities. \emph{Journal of Statistical Theory and
  Practice}, \bold{3}(1), 169--181.
}
\examples{
data("WeatherTask", package = "betareg")
library("flexmix")
wt_betamix <- betamix(agreement ~ 1, data = WeatherTask, k = 2,
  extra_components = extraComponent(type = "betareg", coef =
    list(mean = 0, precision = 2)),
  FLXconcomitant = FLXPmultinom(~ priming + eliciting))
}
\keyword{datasets}
