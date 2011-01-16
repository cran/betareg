\name{betamix}
\encoding{latin1}

\alias{betamix}
\alias{extraComponent}
\alias{fitted,FLXMRbeta-method}
\alias{fitted,betamix-method}
\alias{posterior,betamix,ANY-method}
\alias{clusters,betamix,ANY-method}
\alias{predict,FLXMRbeta-method}
\alias{predict,FLXMRbetafix-method}
\alias{predict,betamix-method}

\title{Finite Mixtures of Beta Regression for Rates and Proportions}

\description{
  Fit finite mixtures of beta regression models for rates and
  proportions via maximum likelihood with the EM algorithm using a
  parametrization with mean (depending through a link function on the
  covariates) and precision parameter (called phi).
}

\usage{
betamix(formula, data, k, fixed, subset, na.action, 
  link = c("logit", "probit", "cloglog", "cauchit", "log",
  "loglog"), link.phi = "log", control = betareg.control(...),
  FLXconcomitant = NULL, extra_components, verbose = FALSE, ID,
  nstart = 1, FLXcontrol = list(), cluster = NULL,
  which = "BIC", ...)

extraComponent(type = c("uniform", "betareg"), coef, delta,
  link = "logit", link.phi = "log")
}

\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}
    or \code{y ~ x | z}; for details see \code{\link{betareg}}).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{k}{a vector of integers indicating the number of components of
    the finite mixture; passed in turn to the \code{k} argument
    of \code{\link[flexmix]{stepFlexmix}}.}
  \item{fixed}{symbolic description of the model for the parameters
    fixed over components (of type \code{~ x | z})}.
  \item{link}{character specification of the link function in
    the mean model (mu). Currently, \code{"logit"}, \code{"probit"},
    \code{"cloglog"}, \code{"cauchit"}, \code{"log"}, \code{"loglog"} are supported.
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{link.phi}{character specification of the link function in
    the precision model (phi). Currently, \code{"identity"},
    \code{"log"}, \code{"sqrt"} are supported. The default is \code{"log"}
    unless \code{formula} is of type \code{y ~ x} where the default is
    \code{"identity"} (for backward compatibility).
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{control}{a list of control arguments specified via
    \code{\link{betareg.control}}.}
  \item{FLXconcomitant}{concomitant variable model; object of class
    \code{FLXP}. Default is the object returned by calling
    \code{\link[flexmix]{FLXPconstant}}.}
  \item{extra_components}{a list containing objects returned by
    \code{extraComponent()}.}
  \item{verbose}{a logical; if \code{TRUE} progress information is shown
    for different starts of the EM algorithm.}
  \item{ID}{grouping variable indicating if observations are from the same
    individual, i.e. the component membership is restricted to be the
    same for these observations.}
  \item{nstart}{for each value of \code{k} run
    \code{\link[flexmix]{stepFlexmix}} \code{nstart} times and keep only
    the solution with maximum likelihood.}
  \item{FLXcontrol}{object of class \code{"FLXcontrol"} or a named list;
    controls the EM algorithm and passed in turn to the \code{control}
    argument of \code{\link[flexmix]{flexmix}}.}
  \item{cluster}{Either a matrix with \code{k} columns of initial
    cluster membership probabilities for each observation; or a factor
    or integer vector with the initial cluster
    assignments of observations at the start of the EM
    algorithm. Default is random assignment into \code{k}
    clusters.}
  \item{which}{number of model to get if \code{k} is a vector of
    integers longer than one. If character, interpreted as
    number of components or name of an information criterion.}
  \item{\dots}{arguments passed to \code{\link{betareg.control}}.}

  \item{type}{specifies if the component follows a uniform distribution
    or a beta regression model.}
  \item{coef}{a vector with the coefficients to determine the midpoint
    of the uniform distribution or names list with the coefficients for
    the mean and precision of the beta regression model.}
  \item{delta}{numeric; half-length of the interval of the uniform
    distribution.}
}

\details{
  The arguments and the model specification are similar to
  \code{\link{betareg}}. Internally \code{\link[flexmix]{stepFlexmix}}
  is called with suitable arguments to fit the finite mixture model with
  the EM algorithm.

  \code{extra_components} is a list where each element corresponds to a
  component where the parameters are fixed a-priori.
}

\value{
  An object of class \code{"flexmix"} containing the best model with
  respect to the log likelihood or the one selected according to
  \code{which} if \code{k} is a vector of integers longer than 1.
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Grün, B., and Leisch, F. (2008). FlexMix Version 2: Finite Mixtures
  with Concomitant Variables and Varying and Constant Parameters.
  \emph{Journal of Statistical Software}, \bold{28}(4), 1--35.
  \url{http://www.jstatsoft.org/v28/i04/}.

Leisch, F. (2004). FlexMix: A General Framework for Finite Mixture
  Models and Latent Class Regression in R.
  \emph{Journal of Statistical Software}, \bold{11}(8), 1--18.
  \url{http://www.jstatsoft.org/v11/i08/}.
}

\author{
  Bettina Grün and Achim Zeileis
}

\seealso{
  \code{\link{betareg}}, \code{\link[flexmix]{flexmix}},
  \code{\link[flexmix]{stepFlexmix}}
}

\examples{
## data with two groups of dyslexic and non-dyslexic children
data("ReadingSkills", package = "betareg")

## try to capture accuracy ~ iq relationship (without using dyslexia
## information) using two beta regression components and one additional
## extra component for a perfect reading score
rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
  nstart = 10, extra_components = extraComponent(type = "uniform",
  coef = 0.99, delta = 0.01))

## visualize result
## intensities based on posterior probabilities
prob <- 2 * (posterior(rs_mix)[cbind(1:nrow(ReadingSkills), clusters(rs_mix))] - 0.5)
## associated HCL colors
col0 <- hcl(c(130, 260, 0), 65, 45, fixup = FALSE)
col1 <- col0[clusters(rs_mix)]
col2 <- hcl(c(130, 260, 0)[clusters(rs_mix)], 65 * abs(prob)^1.5, 95 - 50 * abs(prob)^1.5, fixup = FALSE)
## scatter plot
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19, cex = 1.5, xlim = c(-2, 2))
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
## fitted lines
iq <- -30:30/10
cf <- rbind(c(qlogis(0.99), 0), coef(rs_mix, model = "mean", component = 2:3))
for(i in 1:3) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col0[i]) 
}

\keyword{regression}
\keyword{cluster}
