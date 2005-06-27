\name{envelope.beta}
\alias{envelope.beta}
\title{ Half-Normal Plot of Standardized and Deviance Residuals }
\description{
  Two plots are provided together with an envelope: A Half-Normal plot of standardized residuals and a Half-Normal plot of
  deviance residuals.
}
\usage{
envelope.beta(model=fit.model,sim=100,conf=.90, pch="+",font.main=1, cex.main=1.5, type = c("standardized","deviance"))
}
\arguments{
  \item{model}{ Fitted model using betareg }
  \item{sim}{ The number of simulations used to produce the confidence bounds }
  \item{conf}{ The confidence level }
  \item{pch}{Determines if the points are presented as "+" or "-"}
  \item{font.main}{Determines the type of the font in the main title}
  \item{cex.main}{The size of the font in the main title}
  \item{type}{Determines if the type of the residual is "standardized" or "deviance"}
}
\author{ Alexandre de Bustamante Simas, Department of Statistics, Federal University of Pernambuco (alexandre@cox.de.ufpe.br).

Thanks to Elias Krainski.}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.
}
\seealso{\code{\link{betareg}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}}  }
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit <- betareg(y ~ X)
envelope.beta(fit)
}
\keyword{regression}
