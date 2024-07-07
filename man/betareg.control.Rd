\name{betareg.control}
\alias{betareg.control}

\title{Control Parameters for Beta Regression}
\description{
  Various parameters that control fitting of beta regression models
  using \code{\link{betareg}}.
}

\usage{
betareg.control(phi = TRUE, method = "BFGS", maxit = 5000,
  gradient = NULL, hessian = FALSE, trace = FALSE, start = NULL,
  fsmaxit = 200, fstol = 1e-8, quad = 20, \dots)
}

\arguments{
  \item{phi}{logical indicating whether the precision parameter
    phi should be treated as a full model parameter (\code{TRUE}, default)
    or as a nuisance parameter.}
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}. Additionally, \code{method = "nlminb"}
    can be used to employ \code{\link[stats]{nlminb}}, instead.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{trace}{logical or integer controlling whether tracing information on  
    the progress of the optimization should be produced (passed to \code{\link[stats]{optim}}).}
  \item{gradient}{logical. Should the analytical gradient be used for optimizing
    the log-likelihood? If set to \code{FALSE} a finite-difference approximation
    is used instead. The default of \code{NULL} signals that analytical gradients
    are only used for the classical \code{"beta"} distribution but not for
    \code{"xbetax"} or \code{"xbeta"}.}
  \item{hessian}{logical. Should the numerical Hessian matrix from the \code{optim} output
    be used for estimation of the covariance matrix? By default the analytical solution is employed.
    For details see below.}
  \item{start}{an optional vector with starting values for all parameters (including phi).}
  \item{fsmaxit}{integer specifying maximal number of additional (quasi) Fisher scoring
    iterations. For details see below.}
  \item{fstol}{numeric tolerance for convergence in (quasi) Fisher scoring.
    For details see below.}
  \item{quad}{numeric. The number of quadrature points for numeric
    integration in case of \code{dist = "xbetax"} is used in the beta regression.}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{betareg}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{betareg.control}}.
  Most arguments are passed on directly to \code{optim}, and \code{start} controls
  how \code{optim} is called.

  After the \code{optim} maximization, an additional (quasi) Fisher scoring
  can be perfomed to further enhance the result or to perform additional bias reduction.
  If \code{fsmaxit} is greater than zero, this additional optimization is
  performed and it converges if the threshold \code{fstol} is attained
  for the cross-product of the step size.

  Starting values can be supplied via \code{start} or estimated by
  \code{\link[stats:lmfit]{lm.wfit}}, using the link-transformed response.
  Covariances are in general derived analytically. Only if \code{type = "ML"} and
  \code{hessian = TRUE}, they are determined numerically using the Hessian matrix
  returned by \code{optim}. In the latter case no Fisher scoring iterations are
  performed.
  
  The main parameters of interest are the coefficients in the linear predictor of the
  model and the additional precision parameter phi which can either
  be treated as a full model parameter (default) or as a nuisance parameter. In the latter case
  the estimation does not change, only the reported information in output from \code{print},
  \code{summary}, or \code{coef} (among others) will be different. See also examples.
}

\value{
  A list with the arguments specified.
}

\seealso{\code{\link{betareg}}}

\examples{
options(digits = 4)

data("GasolineYield", package = "betareg")

## regression with phi as full model parameter
gy1 <- betareg(yield ~ batch + temp, data = GasolineYield)
gy1

## regression with phi as nuisance parameter
gy2 <- betareg(yield ~ batch + temp, data = GasolineYield, phi = FALSE)
gy2

## compare reported output
coef(gy1)
coef(gy2)
summary(gy1)
summary(gy2)
}

\keyword{regression}
