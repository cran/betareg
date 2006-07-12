\name{predict.betareg}
\alias{predict.betareg}
\alias{predict}
\title{ Predicted values from beta regression model. }
\description{
  This function returns predictions from a fitted betareg object..}
\usage{
predict.betareg(object, newdata = NULL, type = c("link", "response"), ... )
}
\arguments{
  \item{object}{ fitted model using betareg }
  \item{newdata}{optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}
  \item{type}{the type of prediction required.  The default is on the scale of the linear predictors; the alternative '"response"' is on the scale of the response variable.}
  \item{...}{Optional arguments}
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br). 
This version was modified to fit the standards of the 'lm' and 'glm' models. Andréa Vanessa Rocha (andrea@cox.de.ufpe.br) helped in the development of this new function. }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
summary(fit)
new <- data.frame(x = cbind(0,1,1,0,1,0,1,1,1,201))
predict(fit)
predict(fit,new)
predict(fit,new,type="response")
}
}
\keyword{regression}
