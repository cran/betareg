\name{summary.betareg}
\alias{summary.betareg}
\title{Summary method for Beta Regression}
\description{
  This function returns the coefficient estimates, the precision parameter phi estimate and the Pseudo R2.}
\usage{
summary.betareg(object,...)
}
\arguments{
  \item{object}{ Fitted model using betareg }
  \item{...}{Optional arguments}
}
\value{
  \item{coefficients}{a p by 4 matrix with the columns for the estimated coefficients, the std. errors, the critical values of a normal distribution and the corresponded p-values.}
  \item{estimated phi}{the estimated value of the precision parameter.}
  \item{Pseudo R2}{the value of the pseudo R2 which is described by Ferrari and Cribari-Neto(2004).}
}
\author{ Alexandre de Bustamante Simas, Department of Statistics, Federal University of Pernambuco (alexandredbs@cox.de.ufpe.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, forthcoming.
}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit <- betareg(y ~ X)
summary(fit)
  }
}
\keyword{regression}
