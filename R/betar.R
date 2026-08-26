## beta distribution in regression parameterization (BetaR)
## (mean = mu, precision = phi, support = (0, 1))

dbetar <- function(x, mu, phi, log = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  out <- dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log)
  out[x <= 0 | x >= 1] <- if(log) -Inf else 0
  return(out)
}

pbetar <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  pbeta(q, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
}

qbetar <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
}

rbetar <- function(n, mu, phi) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
}

sbetar <- function(x, mu, phi, which = NULL, drop = TRUE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  p <- c("mu", "phi")
  if (is.null(which)) which <- p
  which <- match.arg(which, p, several.ok = TRUE)
  xstar <- qlogis(x)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ## compute scores
  scr <- function(par) switch(par,
    "mu"  = phi * (xstar - mustar),
    "phi" = mu * (xstar - mustar) + log(1 - x) - digamma((1 - mu) * phi) + digamma(phi))

  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    s <- scr(which)
  } else {
    s <- lapply(which, scr)
    s <- do.call("cbind", s)
    colnames(s) <- which
  }
  return(s)
}

hbetar <- function(x, mu, phi, which = NULL, drop = TRUE, expected = TRUE) {
  if (!expected) NextMethod()

  ## available and selected parameters/combinations and mappings for symmetries
  p <- c("mu" = "mu", "phi:mu" = "mu:phi", "mu:phi" = "mu:phi", "phi" = "phi")
  if (is.null(which)) which <- names(p)
  
  ## which combinations need to be computed?
  which <- match.arg(which, names(p), several.ok = TRUE)
  w <- unique(p[which])

  ## sanity checks
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  n <- max(length(x), length(mu), length(phi))
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)

  ## function for computing Hessian elements (expected only)
  hess <- function(par) switch(par,
    "mu"    = -phi^2 * (psi1 + psi2),
    "phi" = -(psi1 * mu^2 + psi2 * (1 - mu)^2 - trigamma(phi)),
    -phi * (mu * psi1 - (1 - mu) * psi2))
  
  ## if possible return single vector, otherwise collect in matrix
  if (drop && length(which) == 1L) {
    h <- hess(w)
  } else {
    h <- lapply(w, hess)
    h <- do.call("cbind", h)
    colnames(h) <- w
    if (!identical(w, which)) h <- h[, p[which], drop = FALSE]
    colnames(h) <- which
  }
  return(h)
}


## distributions3 interface

BetaR <- function(mu = numeric(), phi = numeric()) {
  n <- c(length(mu), length(phi))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0)
  )
  d <- data.frame(mu = mu, phi = phi)
  class(d) <- c("BetaR", "distribution")
  d
}

mean.BetaR <- function(x, ...) {
  setNames(x$mu, names(x))
}

variance.BetaR <- function(x, ...) {
  rval <- x$mu * (1 - x$mu)/(1 + x$phi)
  setNames(rval, names(x))
}

skewness.BetaR <- function(x, ...) {
  a <- x$mu * x$phi
  b <- (1 - x$mu) * x$phi
  rval <- (6 * ((a - b)^2 * (a + b + 1) - (a * b) * (a + b + 2))) / (a * b * (a + b + 2) * (a + b + 3))
  setNames(rval, names(x))
}

kurtosis.BetaR <- function(x, ...) {
  a <- x$mu * x$phi
  b <- (1 - x$mu) * x$phi
  rval <- (6 * ((a - b)^2 * (a + b + 1) - (a * b) * (a + b + 2))) / (a * b * (a + b + 2) * (a + b + 3))
  setNames(rval, names(x))
}

random.BetaR <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rbetar(n = at, mu = d$mu, phi = d$phi)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.BetaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbetar(x = at, mu = d$mu, phi = d$phi, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.BetaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbetar(x = at, mu = d$mu, phi = d$phi, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.BetaR <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pbetar(q = at, mu = d$mu, phi = d$phi, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.BetaR <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qbetar(p = at, mu = d$mu, phi = d$phi, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.BetaR <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(rep.int(0, length(d)), rep.int(1, length(d)), d, drop = drop)
}

is_discrete.BetaR <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.BetaR <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}

score.BetaR <- function(d, x, which = NULL, drop = TRUE, ...) {
  s <- sbetar(x, mu = d$mu, phi = d$phi, which = which, drop = drop)
  if (!is.null(nam <- names(d))) {
    if (is.null(dim(s))) {
      names(s) <- nam
    } else {
      rownames(s) <- nam
    }
  }
  return(s)
}

hessian.BetaR <- function(d, x, which = NULL, drop = TRUE, expected = TRUE, ...) {
  if (!expected) {
    h <- NextMethod()
  } else {
    h <- hbetar(x, mu = d$mu, phi = d$phi, which = which, drop = drop, expected = expected)
    if (!is.null(nam <- names(d))) {
      if (is.null(dim(h))) {
        names(h) <- nam
      } else {
        rownames(h) <- nam
      }
    }
  }
  return(h)
}
