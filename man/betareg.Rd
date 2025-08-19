\name{betareg}
\encoding{utf8}

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
  link.phi = NULL, type = c("ML", "BC", "BR"), dist = NULL, nu = NULL,
  control = betareg.control(\dots), model = TRUE,
  y = TRUE, x = FALSE, \dots)

betareg.fit(x, y, z = NULL, weights = NULL, offset = NULL,
  link = "logit", link.phi = "log", type = "ML", control = betareg.control(),
  dist = NULL, nu = NULL)
}

\arguments{
  \item{formula}{symbolic description of the model, either of type \code{y ~ x}
    (mean submodel, constant precision) or \code{y ~ x | z} (submodels for both
    mean and precision); for details see below.}
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
  \item{dist}{character specification of the response distribution.
    Usually, this does not have to be set by the user because by default
    the classical \code{"beta"} distribution is used when all observations
    for the dependent variable are in (0, 1). In the presence of boundary
    observations (0 or 1, which cannot be accomodated by \code{"beta"}) the
    extended-support beta mixture distribution (\code{"xbetax"}) is used.
    Additionally, \code{dist = "xbeta"} can be used with fixed exceedence
    parameter \code{nu}, mostly for testing and debugging purposes.}
  \item{nu}{numeric. The fixed value of the expected exceedence parameter \code{nu} 
    in case the extended-support beta mixture distribution is used. By default,
    \code{nu} does not need to be specified and is estimated if needed. So
    setting \code{nu} is mostly for profiling and debugging.}
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
  beta-distributed with parametrization using mean and precision parameter (called mu and
  phi, respectively). The mean mu is linked, as in generalized linear models (GLMs), to the
  explanatory variables through a link function and a linear predictor. Additionally, the
  precision parameter phi can be linked to another (potentially overlapping) set of
  regressors through a second link function, resulting in a model with variable dispersion
  (see Cribari-Neto and Zeileis 2010).
  Estimation is performed by default using maximum likelihood (ML) via \code{\link[stats]{optim}} with
  analytical gradients and starting values from an auxiliary linear regression
  of the transformed response. Subsequently, the \code{optim} result may be enhanced
  by an additional Fisher scoring iteration using analytical gradients and expected information.
  Alternative estimation methods are bias-corrected (BC) or bias-reduced (BR)
  maximum likelihood (see Grün, Kosmidis, and Zeileis 2012). For ML and BC the Fisher
  scoring is just a refinement to move the gradients even closer to zero and can be
  disabled by setting \code{fsmaxit = 0} in the control arguments. For BR the Fisher scoring
  is needed to solve the bias-adjusted estimating equations.

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
  
  Kosmidis and Zeileis (2025) introduce a generalization of the classic beta regression
  model with extended support [0, 1].
  Specifically, the extended-support beta distribution (\code{"xbeta"}) leverages an underlying
  symmetric four-parameter beta distribution with exceedence parameter nu
  to obtain support [-nu, 1 + nu] that is subsequently censored to [0, 1] in order
  to obtain point masses at the boundary values 0 and 1. The extended-support
  beta mixture distribution (\code{"xbetax"}) is a continuous mixture of extended-support
  beta distributions where the exceedence parameter follows an exponential distribution
  with mean nu (rather than a fixed value of nu). The latter \code{"xbetax"}
  specification is used by default in case of boundary observations at 0 and/or 1.
  The \code{"xbeta"} specification with fixed nu is mostly for testing and
  debugging purposes.

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

  The main parameters of interest are the coefficients in the linear predictor of the mean
  model. The additional parameters in the precision model (phi) can either
  be treated as full model parameters (default) or as nuisance parameters. In the latter case
  the estimation does not change, only the reported information in output from \code{print},
  \code{summary}, or \code{coef} (among others) will be different. See also \code{\link{betareg.control}}.

  The implemented algorithms for bias correction/reduction follow Kosmidis and Firth (2010).
  Technical note: In case, either bias correction or reduction is requested,
  the second derivative of the inverse link function is required for \code{link} and
  \code{link.phi}. If the two links are specified by their names (as done by default
  in \code{betareg}), then the \code{"link-glm"} objects are enhanced automatically
  by the required additional \code{d2mu.deta} function. However, if a \code{"link-glm"}
  object is supplied directly by the user, it needs to have the \code{d2mu.deta}
  function or, for backward compatibility, \code{dmu.deta}.

  The original version of the package was written by Alexandre B. Simas and Andrea V. Rocha
  (up to version 1.2). Starting from version 2.0-0 the code was rewritten by Achim Zeileis.
}

\value{
  \code{betareg} returns an object of class \code{"betareg"}, i.e., a list with components as follows.
  For classic beta regressions (\code{dist = "beta"}) several elements are lists with the names \code{"mean"}
  and \code{"precision"} for the information from the respective submodels. For extended-support
  beta regressions (\code{dist = "xbetax"} or \code{"xbeta"}), the corresponding names are
  \code{"mu"} and \code{"phi"} because they are not exactly the mean and precision due to the
  censoring in the response variable.
  
  \code{betareg.fit} returns an unclassed list with components up to \code{converged}.
  
  \item{coefficients}{a list with elements \code{"mean"} (or \code{"mu"}) and \code{"precision"} (or \code{"phi"})
    containing the coefficients from the respective submodels and for extended-support beta regressions
    an additional element \code{"nu"},}
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
  \item{link}{a list with elements \code{"mean"} (or \code{"mu"}) and \code{"precision"} (or \code{"phi"})
    containing the link objects for the respective submodels,}
  \item{converged}{logical indicating successful convergence of \code{optim},}
  \item{call}{the original function call,}
  \item{formula}{the original formula,}  
  \item{terms}{a list with elements \code{"mean"} (or \code{"mu"}), \code{"precision"} (or \code{"phi"}) and
    \code{"full"} containing the terms objects for the respective models,}
  \item{levels}{a list with elements \code{"mean"} (or \code{"mu"}), \code{"precision"} (or \code{"phi"}) and
    \code{"full"} containing the levels of the categorical regressors,}
  \item{contrasts}{a list with elements \code{"mean"} (or \code{"mu"}) and \code{"precision"} (or \code{"phi"})
    containing the contrasts corresponding to \code{levels} from the
    respective models,}
  \item{model}{the full model frame (if \code{model = TRUE}),}
  \item{y}{the response proportion vector (if \code{y = TRUE}),}
  \item{x}{a list with elements \code{"mean"} (or \code{"mu"}) and \code{"precision"} (or \code{"phi"})
    containing the model matrices from the respective models
    (if \code{x = TRUE}).}
}

\references{
Cribari-Neto F, Zeileis A (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \doi{10.18637/jss.v034.i02}

Ferrari SLP, Cribari-Neto F (2004).
  Beta Regression for Modeling Rates and Proportions.
  \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.

Grün B, Kosmidis I, Zeileis A (2012).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  \emph{Journal of Statistical Software}, \bold{48}(11), 1--25.
  \doi{10.18637/jss.v048.i11}

Kosmidis I, Firth D (2010).
  A Generic Algorithm for Reducing Bias in Parametric Estimation.
  \emph{Electronic Journal of Statistics}, \bold{4}, 1097--1112.

Kosmidis I, Zeileis A (2025).
  Extended-Support Beta Regression for [0, 1] Responses.
  \emph{Journal of the Royal Statistical Society C}, forthcoming.
  \doi{10.1093/jrsssc/qlaf039}

Simas AB, Barreto-Souza W, Rocha AV (2010).
  Improved Estimators for a General Class of Beta Regression Models.
  \emph{Computational Statistics & Data Analysis}, \bold{54}(2), 348--366.
}

\seealso{\code{\link{summary.betareg}}, \code{\link{predict.betareg}}, \code{\link{residuals.betareg}},
  \code{\link[Formula]{Formula}}}

\examples{
options(digits = 4)

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
