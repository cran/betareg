\name{residuals.betareg}
\alias{residuals.betareg}
\alias{residuals}
\title{ Residuals function for beta regression models. }
\description{
  This function returns the standardized residuals from beta regression models, deviance residuals or the usual residuals. 
}
\usage{
residuals.betareg(object, type=c("standardized", "usual", "deviance",...)
}
\arguments{
  \item{object}{ Fitted model using betareg. }
  \item{type}{ The desired type of residuals. This function returns by default the standardized residuals, also returns the deviance residuals and may return the usual residuals (which is the value minus estimated mean). }
  \item{...}{Optional arguments}
}
\author{ Alexandre de Bustamante Simas, Department of Statistics, Federal University of Pernambuco (alexandredbs@cox.de.ufpe.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, forthcoming.
}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}} }
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit <- betareg(y ~ X)
residuals(fit)
residuals(fit,type="dev")
  }
}
\keyword{regression}
