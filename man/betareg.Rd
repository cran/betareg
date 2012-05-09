\name{betareg}
\encoding{latin1}

\alias{betareg}
\alias{betareg.fit}

\title{Beta Regression for Rates and Proportions}

\description{
  Fit beta regression models for rates and proportions via maximum likelihood
  using a parametrization with mean (depending through a link function on the
  covariates) and precision parameter (called phi).
}

\usage{
betareg(formula, data, subset, na.action, weights, offset,
  link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
  link.phi = NULL, type = c("ML", "BC", "BR"),
  control = betareg.control(\dots), model = TRUE,
  y = TRUE, x = FALSE, \dots)

betareg.fit(x, y, z = NULL, weights = NULL, offset = NULL,
  link = "logit", link.phi = "log", type = "ML", control = betareg.control())
}

\arguments{
  \item{formula}{symbolic description of the model (of type \code{y ~ x}
    or \code{y ~ x | z}; for details see below).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of case weights.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor for the mean. In \code{betareg.fit},
    \code{offset} may also be a list of two offsets for the mean and precision
    equation, respectively.}
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
  \item{type}{character specification of the type of estimator. Currently,
    maximum likelihood (\code{"ML"}), ML with bias correction (\code{"BC"}),
    and ML with bias reduction (\code{"BR"}) are supported.}
  \item{control}{a list of control arguments specified via
    \code{\link{betareg.control}}.}
  \item{model, y, x}{logicals. If \code{TRUE} the corresponding components
    of the fit (model frame, response, model matrix) are returned.
    For \code{betareg.fit}, \code{x} should be a numeric regressor matrix
    and \code{y} should be the numeric response vector (with values in (0,1)).}
  \item{z}{numeric matrix. Regressor matrix for the precision model, defaulting
    to an intercept only.}
  \item{\dots}{arguments passed to \code{\link{betareg.control}}.}
}

\details{
  Beta regression as suggested by Ferrari and Cribari-Neto (2004) and extended
  by Simas, Barreto-Souza, and Rocha (2010) is implemented in \code{betareg}.
  It is useful in situations where the dependent variable is continuous and restricted to
  the unit interval (0, 1), e.g., resulting from rates or proportions. It is modeled to be
  beta-distributed with parametrization using mean and precision parameter (called phi).
  The mean is linked, as in generalized linear models (GLMs), to the responses through a link
  function and a linear predictor. Additionally, the precision parameter phi can be linked
  to another (potentially overlapping) set of regressors through a second link function,
  resulting in a model with variable dispersion.  
  Estimation is performed by maximum likelihood (ML) via \code{\link[stats]{optim}} using
  analytical gradients and (by default) starting values from an auxiliary linear regression
  of the transformed response. Subsequently, the \code{optim} result may be enhanced
  by an additional Fisher scoring iteration using analytical gradients and expected information.
  This slightly improves the optimization by moving the gradients even closer to zero
  (for \code{type = "ML"} and \code{"BC"}) or solving the bias-adjusted estimating equations
  (for \code{type = "BR"}). For the former two estimators, the optional Fisher scoring
  can be disabled by setting \code{fsmaxit = 0} in the control arguments. See
  Cribari-Neto and Zeileis (2010) and Grün et al. (2012) for details.

  In the beta regression as introduced by Ferrari and Cribari-Neto (2004), the mean of
  the response is linked to a linear predictor described by \code{y ~ x1 + x2} using
  a \code{link} function while the precision parameter phi is assumed to be
  constant. Simas et al. (2009) suggest to extend this model by linking phi to an
  additional set of regressors (\code{z1 + z2}, say): In \code{betareg} this can be
  specified in a formula of type \code{y ~ x1 + x2 | z1 + z2} where the regressors
  in the two parts can be overlapping. In the precision model (for phi), the link
  function \code{link.phi} is used. The default is a \code{"log"} link unless no
  precision model is specified. In the latter case (i.e., when the formula is of type
  \code{y ~ x1 + x2}), the \code{"identity"} link is used by default for backward
  compatibility.
  
  Simas et al. (2009) also suggest further extensions (non-linear specificiations,
  bias correction) which are not yet implemented in \code{betareg}. However,
  Kosmidis and Firth (2010) discuss general algorithms for bias correction/reduction,
  both of which are available in \code{betareg} by setting the \code{type} argument
  accordingly. (Technical note: In case, either bias correction or reduction is requested,
  the second derivative of the inverse link function is required for \code{link} and
  \code{link.phi}. If the two links are specified by their names (as done by default
  in \code{betareg}), then the \code{"link-glm"} objects are enhanced automatically
  by the required additional \code{dmu.deta} function. However, if a \code{"link-glm"}
  object is supplied directly by the user, it needs to have the \code{dmu.deta}
  function.)

  The main parameters of interest are the coefficients in the linear predictor of the mean
  model. The additional parameters in the precision model (phi) can either
  be treated as full model parameters (default) or as nuisance parameters. In the latter case
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

  The original version of the package was written by Alexandre B. Simas and Andrea V. Rocha
  (up to version 1.2). Starting from version 2.0-0 the code was rewritten by Achim Zeileis.
}

\value{
  \code{betareg} returns an object of class \code{"betareg"}, i.e., a list with components as follows.
  \code{betareg.fit} returns an unclassed list with components up to \code{converged}.
  \item{coefficients}{a list with elements \code{"mean"} and \code{"precision"}
    containing the coefficients from the respective models,}
  \item{residuals}{a vector of raw residuals (observed - fitted),}
  \item{fitted.values}{a vector of fitted means,}
  \item{optim}{output from the \code{optim} call for maximizing the log-likelihood(s),}
  \item{method}{the method argument passed to the \code{optim} call,}
  \item{control}{the control arguments passed to the \code{optim} call,}
  \item{start}{the starting values for the parameters passed to the \code{optim} call,}
  \item{weights}{the weights used (if any),}
  \item{offset}{a list of offset vectors used (if any),}
  \item{n}{number of observations,}
  \item{nobs}{number of observations with non-zero weights,}
  \item{df.null}{residual degrees of freedom in the null model (constant mean and dispersion),
    i.e., \code{n - 2},}
  \item{df.residual}{residual degrees of freedom in the fitted model,}
  \item{phi}{logical indicating whether the precision (phi) coefficients will be
    treated as full model parameters or nuisance parameters in subsequent calls to
    \code{print}, \code{summary}, \code{coef} etc.,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{vcov}{covariance matrix of all parameters in the model,}
  \item{pseudo.r.squared}{pseudo R-squared value (squared correlation of linear predictor
    and link-transformed response),}
  \item{link}{a list with elements \code{"mean"} and \code{"precision"}
    containing the link objects for the respective models,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{a list with elements \code{"mean"}, \code{"precision"} and
    \code{"full"} containing the terms objects for the respective models,}
  \item{levels}{a list with elements \code{"mean"}, \code{"precision"} and
    \code{"full"} containing the levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"mean"} and \code{"precision"}
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response proportion vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"mean"} and \code{"precision"}
    containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Ferrari, S.L.P., and Cribari-Neto, F. (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Grün, B., Kosmidis, I., and Zeileis, A. (2012).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  \emph{Journal of Statistical Software}, \bold{48}(11), 1--25.
  \url{http://www.jstatsoft.org/v48/i11/}.

Kosmidis, I., and Firth, D. (2010).
  A Generic Algorithm for Reducing Bias in Parametric Estimation.
  \emph{Electronic Journal of Statistics}, \bold{4}, 1097--1112.

Simas, A.B., Barreto-Souza, W., and Rocha, A.V. (2010).
  Improved Estimators for a General Class of Beta Regression Models.
  \emph{Computational Statistics & Data Analysis}, \bold{54}(2), 348--366.
}

\seealso{\code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}},
  \code{\link[Formula]{Formula}}}

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


## Section 3 from online supplements to Simas et al. (2010)
## mean model as in gy above
## precision model with regressor temp
gy2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)

## MLE column in Table 19
summary(gy2)

## LRT row in Table 18
lrtest(gy, gy2)
}

\keyword{regression}
