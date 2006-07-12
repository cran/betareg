\name{df.residual.betareg}
\alias{df.residual.betareg}
\title{Residual Degrees-of-Freedom}
\description{
  Returns the residual degrees-of-freedom extracted from a fitted betareg model object.
  }
\usage{
df.residual.betareg(object, ...)
}
\arguments{
  \item{object}{ an object for which the degrees-of-freedom are desired. }
  \item{...}{ further arguments passed to or from other methods. }
}
\value{
The value of the residual degrees-of-freedom extracted from the
     object 'x'.
}
\author{ Alexandre de Bustamante Simas, Instituto de Matemática Pura e Aplicada (alesimas@impa.br).
}
\seealso{\code{\link{br.fit}}, \code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
fit <- betareg(oil ~ batch1 + batch2 + batch3 + batch4 + batch5 + 
batch6 + batch7 + batch8 + batch9 + temp, data=pratergrouped)
df.residual(fit)
}
\keyword{regression}
