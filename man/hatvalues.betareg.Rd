\name{hatvalues.betareg}
\alias{hatvalues.betareg}
\title{Hat values}
\description{
   Compute the diagonal elements of the "Hat" matrix of a beta regression model.
  }
\usage{
hatvalues.betareg(model, ...)
}
\arguments{
  \item{model}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
This function returns a vector containing the diagonal elements of the "hat" matrix.
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
hatvalues(fit)
}
\keyword{regression}
