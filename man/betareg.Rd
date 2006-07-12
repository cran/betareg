\name{betareg}
\alias{betareg}
\title{Fitting beta regression models}
\description{
'betareg' is used to fit a regression model where the response is beta distributed using a parameterization of the beta law that is indexed by mean and dispersion parameters. 
  }
\usage{
betareg(formula, link = "logit", data)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit: response ~ covariates. }
  \item{link}{ a link function; the following links are available: logit, probit and cloglog. The default link function is logit. }
  \item{data}{ the 'data argument': a data frame containing the data. }
}
\details{
The model is useful for situations where the variable of interest is continuous, restricted to the standard unit interval (0,1), and related to other variables through a regression structure. The regression parameters of the beta regression model are interpretable in terms of the mean of the response and, when the logit link is used, of an odds ratio, unlike the parameters of a linear regression that employs a transformed response. Estimation is performed by maximum likelihood. The log-likelihood function is maximized using the quasi-Newton BFGS algorithm with analytical first derivatives; the choice of initial values follows the proposal made by Ferrari and Cribari-Neto (2004). 
}
\value{
The function returns an object of class 'betareg' or 'lm'. The function 'summary' is used to obtain an estimate of the precision parameter (phi), and the pseudo R2 value.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit1 <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
fit2 <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + batch6 + 
batch7 + batch8 + batch9 + temp, link = "probit", data=pratergrouped)
summary(fit1)
par(mfrow=c(2,2))
plot(fit1)
}
\keyword{regression}
