\name{GasolineYield}
\alias{GasolineYield}

\title{Estimation of Gasoline Yields from Crude Oil}

\description{
Operational data of the proportion of crude oil converted to gasoline after
distillation and fractionation.
}

\usage{data("GasolineYield")}

\format{
A data frame containing 32 observations on 6 variables.
\describe{
  \item{yield}{proportion of crude oil converted to gasoline after distillation and fractionation.}
  \item{gravity}{crude oil gravity (degrees API).}
  \item{pressure}{vapor pressure of crude oil (lbf/in2).}
  \item{temp10}{temperature (degrees F) at which 10 percent of crude oil has vaporized.}
  \item{temp}{temperature (degrees F) at which all gasoline has vaporized.}
  \item{batch}{factor indicating unique batch of conditions \code{gravity},
    \code{pressure}, and \code{temp10}.}
}
}

\details{
This dataset was collected by Prater (1956), its dependent variable is the
proportion of crude oil after distillation and fractionation. This dataset was
analyzed by Atkinson (1985), who used the linear regression model and noted that
there is ``indication that the error distribution is not quite symmetrical,
giving rise to some unduly large and small residuals'' (p. 60).

The dataset contains 32 observations on the response and on the independent
variables. It has been noted (Daniel and Wood, 1971, Chapter 8) that there are only
ten sets of values of the first three explanatory variables which correspond to
ten different crudes and were subjected to experimentally controlled distillation
conditions. These conditions are captured in variable \code{batch} and
the data were ordered according to the ascending order of \code{temp10}.
}

\source{
Taken from Prater (1956).
}

\references{
Atkinson, A.C. (1985).
  \emph{Plots, Transformations and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis}.
  New York: Oxford University Press.

Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Daniel, C., and Wood, F.S. (1971).
  \emph{Fitting Equations to Data}.
  New York: John Wiley and Sons.

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Prater, N.H. (1956).
  Estimate Gasoline Yields from Crudes.
  \emph{Petroleum Refiner}, \bold{35}(5), 236--238.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy1 <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)
summary(gy1)

## Ferrari and Cribari-Neto (2004)
gy2 <- betareg(yield ~ batch + temp, data = GasolineYield)
## Table 1
summary(gy2)
## Figure 2
par(mfrow = c(3, 2))
plot(gy2, which = 1, type = "pearson", sub.caption = "")
plot(gy2, which = 1, type = "deviance", sub.caption = "")
plot(gy2, which = 5, type = "deviance", sub.caption = "")
plot(gy2, which = 4, type = "pearson", sub.caption = "")
plot(gy2, which = 2:3)
par(mfrow = c(1, 1))

## exclude 4th observation
gy2a <- update(gy2, subset = -4)
gy2a
summary(gy2a)
}

\keyword{datasets}
