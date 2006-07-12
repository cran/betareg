\name{logLik.betareg}
\alias{logLik.betareg}
\title{Extract Log-Likelihood}
\description{
  Obtain the log-likelihood value of a fitted betareg model.
  }
\usage{
logLik.betareg(object, ...)
}
\arguments{
  \item{object}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
Returns the log-likehood value of the model together with the number of parameters.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
logLik(fit)
}
\keyword{regression}
