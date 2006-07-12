\name{cooks.distance.betareg}
\alias{cooks.distance.betareg}
\title{Cook's distance values}
\description{
  Compute the Cook's distance values for beta regression models.
  }
\usage{
cooks.distance.betareg(model, ...)
}
\arguments{
  \item{model}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
This function returns a vector containing the Cook's distance values. 
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.

COOK, R. D. and WEISBERG, S. (1982) \emph{Residuals and Influence in Regression}. London: Chapman and Hall.

}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
cooks.distance(fit)
}
\keyword{regression}
