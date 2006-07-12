\name{residuals.betareg}
\alias{residuals.betareg}
\alias{residuals}
\title{ Residuals function for beta regression models. }
\description{
  This function returns the standardized residuals from beta regression models, deviance residuals or the raw residuals. 
}
\usage{
residuals.betareg(object, type=c("standardized", "raw", "deviance"),...)
}
\arguments{
  \item{object}{ Fitted model using betareg. }
  \item{type}{ The desired type of residuals. This function returns by default the standardized residuals, also returns the deviance residuals and may return the raw residuals (which is the value minus estimated mean). }
  \item{...}{Optional arguments}
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
residuals(fit)
residuals(fit,type="dev")
  }
}
\keyword{regression}
