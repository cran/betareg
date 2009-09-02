\name{betareg}
\alias{betareg}
\alias{betareg.fit}

\title{Beta Regression for Rates and Proportions}

\description{
  Fit beta regression models for rates and proportions via maximum likelihood
  using a parametrization with mean (depending through a link function on the
  covariates) and dispersion parameter (called phi).
}

\usage{
betareg(formula, data, subset, na.action, weights, offset,
  link = c("logit", "probit", "cloglog"),
  control = betareg.control(\dots),
  model = TRUE, y = TRUE, x = FALSE, \dots)

betareg.fit(x, y, weights = NULL, offset = NULL,
  link = "logit", control = betareg.control())
}

\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{link}{character specification of link function.}
  \item{control}{a list of control arguments specified via
    \code{\link{betareg.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.
    For \code{betareg.fit}, \code{x} should be a numeric regressor matrix
    and \code{y} should be the numeric response vector (with values in (0,1)).}
  \item{\dots}{arguments passed to \code{\link{betareg.control}}.}
}

\details{
  Beta regression as suggested by Ferrari and Cribari-Neto (2004) is implemented in \code{betareg}.
  It is useful in situations where the dependent variable is continuous and restricted to
  the unit interval (0, 1), e.g., resulting from rates or proportions. It is modeled to be
  beta-distributed with parametrization using mean and precision/dispersion parameter (called phi).
  The mean is linked, as in generalized linear models (GLMs), to the responses through a link
  function and a linear predictor. Estimation is performed by maximum likelihood (ML) via
  \code{\link[stats]{optim}} using analytical gradients and (by default) starting values from an auxiliary
  linear regression of the transformed response.

  The main parameters of interest are the coefficients
  in the linear predictor and the additional precision/dispersion parameter phi which can either
  be treated as a full model parameter (default) or as a nuisance parameter. In the latter case
  the estimation does not change, only the reported information in output from \code{print},
  \code{summary}, or \code{coef} (among others) will be different. See also \code{\link{betareg.control}}.

  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"betareg"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[graphics]{plot}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{predict}}, \code{\link[stats]{terms}},
  \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}},
  \code{cooks.distance} and \code{hatvalues} (see \code{\link[stats]{influence.measures}}),
  \code{\link[betareg]{gleverage}} (new generic), \code{\link[sandwich]{estfun}} and
  \code{\link[sandwich]{bread}} (from the \pkg{sandwich} package), and
  \code{\link[lmtest]{coeftest}} (from the \pkg{lmtest} package).
  
  See \code{\link{predict.betareg}}, \code{\link{residuals.betareg}}, \code{\link{plot.betareg}},
  and \code{\link{summary.betareg}} for more details on all methods.
}

\value{
  \code{betareg} returns an object of class \code{"betareg"}, i.e., a list with components as follows.
  \code{betareg.fit} returns an unclassed list with components up to \code{converged}.
  \item{coefficients}{vector with estimated regression coefficients and dispersion
    (or precision) parameter phi,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{output from the \code{optim} call for maximizing the log-likelihood(s),}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call,}
  \item{weights}{the weights used (if any),}
  \item{offset}{the offset vector used (if any),}
  \item{n}{number of observations,}
  \item{df.null}{residual degrees of freedom for the null model (= \code{n - 2}),}
  \item{df.residual}{residual degrees of freedom for fitted model,}
  \item{phi}{logical indicating whether phi will be treated as a full model parameter
    or a nuisance parameter in subsequent calls to \code{print}, \code{summary},
    \code{coef} etc.,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of all parameters in the model (including phi),}
  \item{pseudo.R.squared}{pseudo R-squared value (squared correlation of linear predictor
    and link-transformed response),}
  \item{link}{link object used,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{the terms object used,}
  \item{levels}{levels of the categorical regressors,}
  \item{contrasts}{contrasts corresponding to \code{levels},}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response proportion vector (if \code{y = TRUE}),}
  \item{x}{the model matrix (if \code{x = TRUE}).}
}

\references{
Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
}

\seealso{\code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}}}

\examples{
## Section 4 from Ferrari and Cribari-Neto (2004)
data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")

## Table 1
gy <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy)

## Table 2
fe_lin <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)
library("lmtest")
bptest(fe_lin)
fe_beta <- betareg(I(food/income) ~ income + persons, data = FoodExpenditure)
summary(fe_beta)

## nested model comparisons via Wald and LR tests
fe_beta2 <- betareg(I(food/income) ~ income, data = FoodExpenditure)
lrtest(fe_beta, fe_beta2)
waldtest(fe_beta, fe_beta2)
}

\keyword{regression}
