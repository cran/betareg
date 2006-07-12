\name{br.fit}
\alias{br.fit}
\title{ Function to estimate the coefficients of the Beta Regression. }
\description{
  Function to estimate the coefficients of the Beta Regression.
}
\usage{
br.fit(x, y, link)
}
\arguments{
  \item{x}{ the terms matrix }
  \item{y}{ the response matrix }
  \item{link}{ an object returned by the make.link function. }
}
\value{
   The function returns the estimated coefficients of the regression, the estimation of the precision parameter phi, the fitted values, the residuals,the pseudo R2, the standard errors, the critical values of the normal distribution and the p-values of each coefficient and the phi parameter.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.
}
\seealso{\code{\link{betareg}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}}  }
\examples{
x <- cbind(rnorm(10),rnorm(10))
y <- cbind(rbeta(10,1,1))
link = betareg(y ~ x)$funlink
br.fit(x,y,link)
}
\keyword{regression}
