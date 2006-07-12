\name{pratergrouped}
\alias{pratergrouped}
\docType{data}
\title{Grouped Prater's gasoline data}
\description{
  Operational data of the proportion of crude oil converted to gasoline after distilation and fractionation.}

\usage{data(pratergrouped)}
\format{
  A data frame with 32 observations on 11 variables.
  More details on the independent variables can be found in the details section below.

  [,1] 'oil'    The proportion of crude oil converted to gasoline after distilation and fractionation.
  [,2] 'batch1'  The first batch of crude oil.
  [,3] 'batch2'  The second batch of crude oil.
  [,4] 'batch3'  The third batch of crude oil.
  [,5] 'batch4'  The fourth batch of crude oil.
  [,6] 'batch5'  The fifth batch of crude oil.
  [,7] 'batch6'  The sixth batch of crude oil.
  [,8] 'batch7'  The seventh batch of crude oil.
  [,9] 'batch8'  The eighth batch of crude oil.
  [,10] 'batch9' The ninth batch of crude oil.
  [,11] 'temp'   The temperature (degrees F) at which all the gasoline is vaporized.



}
\details{
The dataset contains 32 observations on the response and on the independent variables. It has been noted (Daniel and Wood, 1971, Ch. 8) that there are only ten sets of values of the first three explanatory variables which correspond to ten different crudes and were subjected to experimentally controlled distilation conditions. This dataset was analyzed by Atkinson (1985), who used the linear regression model and noted that there is ``indication that the error distribution is not quite symmetrical, giving rise to some unduly large and small residuals'' (p. 60). He proceeded to transform the response so that the transformed dependent variable assumed values on the real line, and then used it in a linear regression analysis.

The data were ordered according to the ascending order of the covariate that measures the temperature at which 10% of the crude oil has become vapor. This variable assumes ten different values and they are used to define the ten batches of crude oil. The model specification for the mean of the response uses an intercept (x1 = 1), nine dummy variables for the first nine batches of crude oil (x2 , . . . , x10 ) and the covariate that measures the temperature (degrees F) at which all the gasoline is vaporized (x11 ).
}
\source{
  Prater, N.H. (1956). Estimate gasoline yields from crudes. New York: Springer-Verlag.
}
\references{
Atkinson, A.C. (1985). \emph{Plots, Transformations and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis.} New York: Oxford University Press.

Daniel, C. and Wood, F.S. (1971). \emph{Fitting Equations to Data}. New York: Wiley
}
\examples{
data(pratergrouped)
fit1 <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
fit2 <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + batch6 + 
batch7 + batch8 + batch9 + temp, link = "probit", data=pratergrouped)
plot(fit1)
}
\keyword{datasets}
