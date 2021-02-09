\name{plot.betareg}
\alias{plot.betareg}

\title{Diagnostic Plots for betareg Objects}

\description{
  Various types of standard diagnostic plots can be produced, involving various types of
  residuals, influence measures etc.
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
  \item{which}{numeric. If a subset of the plots is required, specify a subset of the numbers \code{1:6}.}
  \item{caption}{character. Captions to appear above the plots.}
  \item{sub.caption}{character. Common title-above figures if there are multiple.}
  \item{main}{character. Title to each plot in addition to the above \code{caption}.}
  \item{ask}{logical. If \code{TRUE}, the user is asked before each plot.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
  \item{type}{character indicating type of residual to be used, see \code{\link{residuals.betareg}}.}
  \item{nsim}{numeric. Number of simulations in half-normal plots.}
  \item{level}{numeric. Confidence level in half-normal plots.}
}

\details{
  The \code{plot} method for \code{\link{betareg}} objects produces various types
  of diagnostic plots. Most of these are standard for regression models and involve
  various types of residuals, influence measures etc. See Ferrari and Cribari-Neto (2004)
  for a discussion of some of these displays.
  
  The \code{which} argument can be used to select a subset of currently six supported
  types of displays. The corresponding element of \code{caption} contains a brief
  description. In some more detail, the displays are: Residuals (as selected by
  \code{type}) vs indices of observations (\code{which = 1}). Cook's distances
  vs indices of observations (\code{which = 2}). Generalized leverage vs
  predicted values (\code{which = 3}). Residuals vs linear predictor  (\code{which = 4}).
  Half-normal plot of residuals (\code{which = 5}), which is obtained using a simulation
  approach. Predicted vs observed values (\code{which = 6}).
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \doi{10.18637/jss.v034.i02}

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
