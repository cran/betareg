\name{pratergrouped}
\alias{pratergrouped}
\docType{data}
\title{Grouped Prater's gasoline data}
\description{
  This dataset was collected by Prater (1956), its dependent variable is the proportion of crude oil after distilation and fractionation. This dataset was analyzed by Atkunson (1985), who used the linear regression model and noted that there is "indication that the error distribution is not quite symmetrical, giving rise to some unduly large and small residuals" (p. 60).  
}
\usage{data(prater)}
\format{
  A data frame with 32 observations on 11 variables.
}
\source{
  Prater, N.H. (1956). Estimate gasoline yields from crudes. New York: Springer-Verlag.
}
\references{
  Atkinson, A.C. (1985). \emph{Plots, Transformations and Regression: An Introduction to Graphical Methods of Diagnostic Regression Analysis.} New York: Oxford University Press.
}
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit1 <- betareg(y ~ X)
fit2 <- betareg(y ~ X, link = "probit")
plot(fit1)
}
\keyword{datasets}
