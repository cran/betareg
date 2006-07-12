\name{gen.lev.betareg}
\alias{gen.lev.betareg}
\title{Generalized leverage values}
\description{
   Compute the generalized leverages values for beta regression models.
  }
\usage{
gen.lev.betareg(x, ...)
}
\arguments{
  \item{x}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
This function returns a vector containing the generalized leverage values.
The generalized leverage was proposed by Wei, Hu and Fung and was adapted to the beta regression model by Ferrari and Cribari-Neto (2004).
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
}
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, v. 31, n. 7, p. 799-815.

WEI, B.-C., HU, Y.-Q. and FUNG, W.-K. (1998). Generalized leverage and its applications. Scandinavian Journal of Statistics, 25, 25-37
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
gen.lev.betareg(fit)
}
\keyword{regression}
