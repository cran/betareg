\name{vcov.betareg}
\alias{vcov.betareg}
\title{Calculate Variance-Covariance Matrix for a betareg Fitted Model Object}
\description{
  Returns the variance-covariance matrix of the main parameters of a betareg fitted model object.
  }
\usage{
vcov.betareg(object, ...)
}
\arguments{
  \item{object}{ a fitted model using betareg. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
A matrix of the estimated covariances between the parameter estimates in the linear or non-linear predictor of the model.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
This function was inspired on the function made by Achim Zeileis.
}
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
vcov(fit)
}
\keyword{regression}
