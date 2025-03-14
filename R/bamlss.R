betar_family <- function(link = "logit", link.phi = "log", ...)
{
  ## handle links as in betareg()
  link <- match.arg(link, c("logit", "probit", "cloglog", "cauchit", "log", "loglog"))
  link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))
  links <- c(mu = link, phi = link.phi)
  
  ## set up family list for bamlss/gamlss2
  rval <- list(
    "family" = "xbetax",
    "names" = names(links),
    "links" =  links,
    "valid.response" = function(x) {
      ok <- all(x > 0 & x < 1)
      if(!ok) stop("response values not in (0, 1)", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) dbetar(y, mu = par$mu, phi = par$phi, log = log),
    "p" = function(y, par, ...) pbetar(y, mu = par$mu, phi = par$phi, ...),
    "q" = function(p, par) qbetar(p, mu = par$mu, phi = par$phi),
    "r" = function(n, par) rbetar(n, mu = par$mu, phi = par$phi),
    "mean" = function(par) par$mu,
    "variance" = function(par) (par$mu * (1 - par$mu))/(1 + par$phi),
    "score" = list(
      "mu"  = function(y, par, ...) sbetar(y, mu = par$mu, phi = par$phi, parameter = "mu",  drop = TRUE),
      "phi" = function(y, par, ...) sbetar(y, mu = par$mu, phi = par$phi, parameter = "phi", drop = TRUE)
    ),
    "hess" = list(
      "mu"  = function(y, par, ...) hbetar(y, mu = par$mu, phi = par$phi, parameter = "mu",  drop = TRUE),
      "phi" = function(y, par, ...) hbetar(y, mu = par$mu, phi = par$phi, parameter = "phi", drop = TRUE)
    )
  )
  class(rval) <- "family.bamlss"
  rval
}

xbetax_family <- function(link = "logit", link.phi = "log", link.nu = "log", quad = 20, tol = .Machine$double.eps^0.7, ...)
{
  ## handle links as in betareg()
  link <- match.arg(link, c("logit", "probit", "cloglog", "cauchit", "log", "loglog"))
  link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))
  link.nu <- match.arg(link.phi, c("identity", "log", "sqrt"))
  links <- c(mu = link, phi = link.phi, nu = link.nu)

  ## set up family list for bamlss/gamlss2
  rval <- list(
    "family" = "xbetax",
    "names" = names(links),
    "links" =  links,
    "valid.response" = function(x) {
      ok <- all(x >= 0 & x <= 1)
      if(!ok) stop("response values not in [0, 1]", call. = FALSE)
      ok
    },
    "d" = function(y, par, log = FALSE) dxbetax(y, mu = par$mu, phi = par$phi, nu = par$nu, log = log, quad = quad),
    "p" = function(y, par, ...) pxbetax(y, mu = par$mu, phi = par$phi, nu = par$nu, quad = quad, ...),
    "q" = function(p, par) qxbetax(p, mu = par$mu, phi = par$phi, nu = par$nu, quad = quad, tol = tol),
    "r" = function(n, par) rxbetax(n, mu = par$mu, phi = par$phi, nu = par$nu),
    "mean" = function(par) mean_xbetax(mu = par$mu, phi = par$phi, nu = par$nu, quad = quad),
    "variance" = function(par) var_xbetax(mu = par$mu, phi = par$phi, nu = par$nu, quad = quad)
  )
  class(rval) <- "family.bamlss"
  rval
}
