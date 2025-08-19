\name{betar_family}
\alias{betar_family}
\alias{xbetax_family}

\title{Family Objects for (Extended-Support) Beta Regression}

\description{
Family objects for distributional regression with the (extended-support)
beta distribution via \pkg{bamlss} or \pkg{gamlss2}.
}

\usage{
betar_family(link = "logit", link.phi = "log", \dots)

xbetax_family(link = "logit", link.phi = "log", link.nu = "log",
  quad = 20, tol = .Machine$double.eps^0.7, \dots)
}

\arguments{
  \item{link}{character specification of the link function in
    the mean model (mu). Currently, \code{"logit"}, \code{"probit"},
    \code{"cloglog"}, \code{"cauchit"}, \code{"log"}, \code{"loglog"} are supported.}
  \item{link.phi}{character specification of the link function in
    the precision model (phi). Currently, \code{"log"}, \code{"identity"},
    \code{"sqrt"} are supported.}
  \item{link.nu}{character specification of the link function in
    the exceedence model (nu). Currently, \code{"log"}, \code{"identity"},
    \code{"sqrt"} are supported.}
  \item{quad}{numeric. The number of quadrature points for numeric
    integration of the continuous mixture in \code{\link{dxbetax}}. Alternatively,
    a matrix with nodes and weights for the quadrature points can be specified.}
  \item{tol}{numeric. Accuracy (convergence tolerance) for numerically
    determining quantiles based on \code{\link[stats]{uniroot}} and
    \code{\link{pxbetax}}.}
  \item{\dots}{Arguments passed to functions that are called within the family object.}
}

\details{
Family objects for \pkg{bamlss} (Umlauf et al. 2019, 2021) and \pkg{gamlss2}
(Umlauf et al. 2025) are essentially lists of functions providing a standardized
interface to the d/p/q/r functions of distributions. Hence, \code{betar_family}
interfaces the classical beta distribution in regression specification, see
\code{\link{dbetar}}. Analogously, \code{xbetax_family} interfaces the
extended-support beta mixture specification (Kosmidis and Zeileis 2025),
see \code{\link{dxbetax}}.
}

\value{
A list of class \code{\link[bamlss]{family.bamlss}}.
}

\seealso{\code{\link{dbetar}}, \code{\link{dxbetax}}, \code{\link[bamlss]{family.bamlss}}}

\references{
Kosmidis I, Zeileis A (2025).
  Extended-Support Beta Regression for [0, 1] Responses.
  \emph{Journal of the Royal Statistical Society C}, forthcoming.
  \doi{10.1093/jrsssc/qlaf039}

Umlauf N, Klein N, Zeileis A (2019).
  BAMLSS: Bayesian Additive Models for Location, Scale and Shape (and Beyond).
  \emph{Journal of Computational and Graphical Statistics}, \bold{27}(3), 612--627.
  \doi{10.1080/10618600.2017.1407325}

Umlauf N, Klein N, Simon T, Zeileis A (2021).
  bamlss: A Lego Toolbox for Flexible Bayesian Regression (and Beyond).
  \emph{Journal of Statistical Software}, \bold{100}(4), 1--53.
  \doi{10.18637/jss.v100.i04}
  
Umlauf N, Stasinopoulos M, Rigby R, Stauffer R (2025).
  gamlss2: Infrastructure for Flexible Distributional Regression.
  R package version 0.1-0.
  \url{https://gamlss-dev.github.io/gamlss2/}
}

\examples{
## IGNORE_RDIFF_BEGIN
\dontshow{ if(!requireNamespace("bamlss")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("not all packages required for the example are installed")
  } else q() }
}
## package and data
library("betareg")
library("bamlss")
data("ReadingSkills", package = "betareg")

## classical beta regression via ML
rs1 <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq, data = ReadingSkills)

## Bayesian additive model (with low number of iterations to speed up the example)
set.seed(0)
rs2 <- bamlss(accuracy ~ s(iq, by = dyslexia) | dyslexia + iq, data = ReadingSkills,
  family = betar_family(), eps = 1e-7, n.iter = 400, burnin = 100)

## Bayesian model shrinks the effects compared to ML
plot(accuracy ~ iq, data = ReadingSkills, pch = 19, col = dyslexia)
nd <- data.frame(
  iq = rep(-19:20/10, 2),
  dyslexia = factor(rep(c("no", "yes"), each = 40), levels = c("no", "yes"))
)
nd$betareg <- predict(rs1, newdata = nd, type = "response")
nd$bamlss  <- predict(rs2, newdata = nd, type = "parameter", model = "mu")
lines(betareg ~ iq, data = nd, subset = dyslexia == "no",  col = 1, lwd = 2, lty = 1)
lines(betareg ~ iq, data = nd, subset = dyslexia == "yes", col = 2, lwd = 2, lty = 1)
lines(bamlss  ~ iq, data = nd, subset = dyslexia == "no",  col = 1, lwd = 2, lty = 2)
lines(bamlss  ~ iq, data = nd, subset = dyslexia == "yes", col = 2, lwd = 2, lty = 2)
legend("topleft", c("Dyslexia: no", "Dyslexia: yes", "betareg", "bamlss"),
  lty = c(0, 0, 1, 2), pch = c(19, 19, NA, NA), col = c(1, 2, 1, 1), bty = "n")
## IGNORE_RDIFF_END

## xbetax_family(): requires more time due to Gaussian quadrature
## for gamlss2: install.packages("gamlss2", repos = "https://gamlss-dev.R-universe.dev")
}
