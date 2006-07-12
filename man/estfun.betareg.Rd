\name{estfun.betareg}
\alias{estfun.betareg}
\title{Extract Empirical Estimating Functions}
\description{
  Function for extracting the empirical estimating functions of a betareg fitted model.
  }
\usage{
estfun.betareg(x, ...)
}
\arguments{
  \item{x}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{ 
A matrix containing the empirical estimating functions. Typically, this should be an n x k matrix corresponding to n observations and k parameters. The columns should be named as in coef or terms, respectively.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
This function was made by Achim Zeileis.
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.

ZEILEIS, A. (2006). Object-oriented Computation of Sandwich Estimators. Package vignette.
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
estfun.betareg(fit)
}
\keyword{regression}
