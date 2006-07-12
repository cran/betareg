\name{anova.betareg}
\alias{anova.betareg}
\title{A preliminary version for the analisys of variance table}
\description{
  Compute an analysis of variance table for one or two model fits.
  }
\usage{
anova.betareg(object, object2, ...)
}
\arguments{
  \item{object}{ a fitted model using betareg, if there are two models, then this model should be the restricted one }
  \item{object2}{ (optional) a fitted model using betareg}
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
If object2 is missing, an anova table is produced considering the null that the model contain only the intercept. If object2 isn't missing (i.e., a restricted model was placed), 
an anova table is produced considering the null that the second model is true. This "anova table" is constructed using the log-likelihood ratio test suggest by Ferrari and Cribari-Neto (2004).
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
Andréa Vanessa Rocha (andrea@cox.de.ufpe.br) contributed in the development of this function.
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit1 <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
fit2 <- betareg(oil ~ batch1 + batch2 + temp, data=pratergrouped)
## With only one model
anova(fit1)
## With two models, note that the first model is the restricted one.
anova(fit2,fit1)
}
\keyword{regression}
