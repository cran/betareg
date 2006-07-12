\name{prater}
\alias{prater}
\docType{data}
\title{Prater's gasoline data}
\description{
Operational data of the proportion of crude oil converted to gasoline after distilation and fractionation.
}

\usage{data(prater)}
\format{
  A data frame with 32 observations on 5 variables.

  [,1] 'oil'        The proportion of crude oil converted to gasoline after distilation and fractionation.
  [,2] 'apigrav'    The crude oil gravity (degrees API).
  [,3] 'vaporpress' The vapor pressure of the crude oil (lbf/in2 ).
  [,4] 'tempten'    The crude oil 10% point ASTM (i.e., the temperature at which 10% of the crude oil has become vapor).
  [,5] 'temp'       The temperature (degrees F) at which all the gasoline is vaporized.
}
\details{
This dataset was collected by Prater (1956), its dependent variable is the proportion of crude oil after distilation and fractionation. This dataset was analyzed by Atkunson (1985), who used the linear regression model and noted that there is ``indication that the error distribution is not quite symmetrical, giving rise to some unduly large and small residuals'' (p. 60).
}

\source{
Prater, N.H. (1956). Estimate gasoline yields from crudes. New York: Springer-Verlag.
}
\references{
Atkinson, A.C. (1985). \emph{Plots, Transformations and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis.} New York: Oxford University Press.
}
\examples{
data(prater)
fit <- betareg(oil ~ apigrav + vaporpress + tempten + temp, data=prater)
plot(fit)
}
\keyword{datasets}
