\name{plot.betareg}
\alias{plot.betareg}
\alias{plot}
\title{ Plot Diagnostics for an betareg Object }
\description{
  This function returns four plots: a plot of residuals against fitted values, a plot of standardized residuals against fitted values, a generalized leverage plot against fitted values and a plot of Cook's distances versus row labels.}
\usage{
plot.betareg(x, which = 1:4, caption = c("Deviance Residuals vs indices of obs.", 
    "Standardized Residuas vs indices of obs.", "Generalized Leverage vs. Predicted Values", "Cook's distance plot"), 
    panel = points, sub.caption = deparse(x$call), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75)
}
\arguments{
  \item{x}{ Fitted model by betareg. }
  \item{which}{If a subset of the plots is required, specify a subset of the numbers '1:4'.}
  \item{caption}{Captions to appear above the plots.}
  \item{panel}{Panel function.  A useful alternative to 'points' is 'panel.smooth'.}
  \item{sub.caption}{common title-above figures if there are multiple; used as 'sub' (s.'title') otherwise.}
  \item{main}{title to each plot-in addition to the above 'caption'.}
  \item{ask}{logical; if `TRUE', the user is \emph{asked} before each plot, see `par(ask=.)'.}
  \item{...}{optional arguments.}
  \item{id.n}{number of points to be labelled in each plot, starting with the most extreme.}
  \item{labels.id}{vector of labels, from which the labels for extreme points will be chosen.  'NULL' uses observation numbers.}
  \item{cex.id}{magnification of point labels.}
}
\author{ Alexandre de Bustamante Simas, Department of Statistics, Federal University of Pernambuco (alexandredbs@cox.de.ufpe.br). }
\references{
FERRARI, S.L.P., CRIBARI-NETO, F. (2004). Beta regression for modeling rates and proportions. Journal of Applied Statistics, forthcoming.
}
\seealso{\code{\link{betareg}}, \code{\link{br.fit}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}} }
\examples{
data(pratergrouped)
attach(pratergrouped)
y <- V11/100
X <- cbind(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10)
fit <- betareg(y ~ X)
par(mfrow=c(2,2))
plot(fit)
  }
}
\keyword{regression}
