\name{predict.betareg}
\alias{predict.betareg}
\alias{predict}
\title{ Predicted values from beta regression model. }
\description{
  This function returns the predicted values based on beta regression models.}
\usage{
predict.betareg(object, terms,...)
}
\arguments{
  \item{object}{ Fitted model using betareg }
  \item{terms}{ A matrix or a vector containing the values that you want to have the predict. }
  \item{...}{Optional arguments}
}
\author{ Alexandre de Bustamante Simas, Department of Statistics, Federal University of Pernambuco (alexandredbs@cox.de.ufpe.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, forthcoming.
}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit <- betareg(y ~ X)
predict(fit)
predict(fit,c(2,3,4,5,6,7,8,9,10,11))
  }
}
\keyword{regression}
