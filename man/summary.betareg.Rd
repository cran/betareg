\name{summary.betareg}
\alias{summary.betareg}
\title{Summary method for Beta Regression}
\description{
  'summary' method for class 'betareg'.}
\usage{
summary.betareg(object,...)
}
\arguments{
  \item{object}{ Fitted model using betareg }
  \item{...}{Optional arguments}
}
\value{
  \item{call}{the component from object}
  \item{deviance residuals}{a summary of the deviance residuals, see residuals.betareg}
  \item{coefficients}{a p by 4 matrix with the columns for the estimated coefficients, the std. errors, the critical values of a normal distribution and the corresponded p-values.}
  \item{estimated phi}{the estimated value of the precision parameter together with its standard error.}
  \item{null deviance}{the deviance from the model only containing the intercept (but, assuming the estimated phi from the current model)}
  \item{residual deviance}{the deviance from the current model}
  \item{log-likelihood ratio statistic}{the log-likelihood ratio statistic from the full model against the model containing only the intercept}
  \item{Pseudo R2}{the value of the pseudo R2 which is described by Ferrari and Cribari-Neto(2004).}
  \item{AIC}{the corresponding AIC} 
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br). 
This version was modified to fit the standards of the 'lm' and 'glm' models. Andréa Vanessa Rocha (andrea@cox.de.ufpe.br) helped in the development of this new function.}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
summary(fit)
  }
}
\keyword{regression}
