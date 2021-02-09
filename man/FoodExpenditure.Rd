\name{FoodExpenditure}
\alias{FoodExpenditure}

\title{Proportion of Household Income Spent on Food}

\description{
  Data on proportion of income spent on food for a random sample of 38 households in a large US city.
}

\usage{data("FoodExpenditure")}

\format{
A data frame containing 38 observations on 3 variables.
\describe{
  \item{food}{household expenditures for food.}
  \item{income}{household income.}
  \item{persons}{number of persons living in household.}
}
}

\source{
Taken from Griffiths et al. (1993, Table 15.4).
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \doi{10.18637/jss.v034.i02}

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Griffiths, W.E., Hill, R.C., and Judge, G.G. (1993).
  \emph{Learning and Practicing Econometrics}
  New York: John Wiley and Sons.
}

\seealso{\code{\link{betareg}}}

\examples{
data("FoodExpenditure", package = "betareg")

## Ferrari and Cribari-Neto (2004)
## Section 4
fe_lin <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)
library("lmtest")
bptest(fe_lin)

## Table 2
fe_beta <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure)
summary(fe_beta)
}

\keyword{datasets}
