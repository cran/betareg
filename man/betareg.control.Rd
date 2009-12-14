\name{betareg.control}
\alias{betareg.control}

\title{Control Parameters for Beta Regression}
\description{
  Various parameters that control fitting of beta regression models
  using \code{\link{betareg}}.
}

\usage{
betareg.control(phi = TRUE, method = "BFGS", maxit = 5000,
  hessian = FALSE, trace = FALSE, start = NULL, \dots)
}

\arguments{
  \item{phi}{logical indicating whether the precision parameter
    phi should be treated as a full model parameter (\code{TRUE}, default)
    or as a nuisance parameter.}
  \item{method}{characters string specifying the \code{method} argument
    passed to \code{\link[stats]{optim}}.}
  \item{maxit}{integer specifying the \code{maxit} argument (maximal number
    of iterations) passed to \code{\link[stats]{optim}}.}
  \item{trace}{logical or integer controlling whether tracing information on  
    the progress of the optimization should be produced (passed to \code{\link[stats]{optim}}).}
  \item{hessian}{logical. Should the numerical Hessian matrix from the \code{optim} output
    be used for estimation of the covariance matrix? If \code{FALSE} (the default),
    the analytical solution is employed.}
  \item{start}{an optional vector with starting values for all parameters (including phi).}
  \item{\dots}{arguments passed to \code{\link[stats]{optim}}.}
}

\details{
  All parameters in \code{\link{betareg}} are estimated by maximum likelihood
  using \code{\link[stats]{optim}} with control options set in \code{\link{betareg.control}}.
  Most arguments are passed on directly to \code{optim}, only \code{start} controls
  how \code{optim} is called.

  Starting values can be supplied via \code{start} or estimated by
  \code{\link[stats:lmfit]{lm.wfit}}, using the link-transformed response.
  Covariances are derived analytically (if \code{hessian = FALSE}, the default) or
  numerically using the Hessian matrix returned by \code{\link[stats]{optim}}.
  
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
