\name{plot.betareg}
\alias{plot.betareg}

\title{Diagnostic Plots for betareg Objects}

\description{
  Various types of standard diagnostic plots can be produced, involving various types of
  residuals, influence measures etc. Half-normal plots can be produced using a simulation approach.
}

\usage{
\method{plot}{betareg}(x, which = 1:4,
  caption = c("Residuals vs indices of obs.", "Cook's distance plot",
    "Generalized leverage vs predicted values", "Residuals vs linear predictor", 
    "Half-normal plot of residuals", "Predicted vs observed values"),
    sub.caption = paste(deparse(x$call), collapse = "\n"), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    \dots, type = "sweighted2", nsim = 100, level = 0.9)
}

\arguments{
  \item{x}{fitted model object of class \code{"betareg"}.}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers \code{1:8}.}
  \item{caption}{captions to appear above the plots.}
  \item{sub.caption}{common title-above figures if there are multiple.}
  \item{main}{title to each plot-in addition to the above \code{caption}.}
  \item{ask}{logical. If \code{TRUE}, the user is asked before each plot.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
  \item{type}{character indicating type of residual to be used, see \code{\link{residuals.betareg}}.}
  \item{nsim}{number of simulations in half-normal plots.}
  \item{level}{confidence level in half-normal plots.}
}

\references{
Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
}

\seealso{\code{\link{betareg}}}

\examples{
data("GasolineYield", package = "betareg")

gy <- betareg(yield ~ gravity + pressure + temp10 + temp, data = GasolineYield)

par(mfrow = c(3, 2))
plot(gy, which = 1:6)
par(mfrow = c(1, 1))
}

\keyword{regression}
