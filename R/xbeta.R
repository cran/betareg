## extended-domain beta distribution (XBeta)
## based on the censored symmetric four-parameter beta distribution in regression parameterization
## (mean = mu, precision = phi, latent support = (-nu, 1 + nu) censored to [0, 1])

dxbeta <- function(x, mu, phi, nu = 0, log = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )

  ## essentially rely on rescaling as in symmetric four-parameter beta distribution
  out <- dbeta((x + nu) / (1 + 2 * nu), shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log)
  out <- if(log) out - log(1 + 2 * nu) else out/(1 + 2 * nu)

  ## unify lengths of all variables
  n <- length(out)
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## boundary cases
  out[x <= 0] <- pbeta((0 + nu[x <= 0]) / (1 + 2 * nu[x <= 0]), shape1 = (mu * phi)[x <= 0], shape2 = ((1 - mu) * phi)[x <= 0], log.p = log, lower.tail = TRUE)
  out[x >= 1] <- pbeta((1 + nu[x >= 1]) / (1 + 2 * nu[x >= 1]), shape1 = (mu * phi)[x >= 1], shape2 = ((1 - mu) * phi)[x >= 1], log.p = log, lower.tail = FALSE)
  out[x < 0 | x > 1] <- if(log) -Inf else 0

  return(out)
}

pxbeta <- function(q, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )

  ## essentially rely on rescaling as in symmetric four-parameter beta distribution
  out <- pbeta((q + nu) / (1 + 2 * nu), shape1 = mu * phi, shape2 = (1 - mu) * phi, lower.tail = lower.tail, log.p = log.p)

  ## unify lengths of all variables
  n <- length(out)
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## boundary cases
  if(lower.tail) {
    out[q <= 0] <- dxbeta(0, mu = mu[q <= 0], phi = phi[q <= 0], nu = nu[q <= 0], log = log.p)
    out[q <  0] <- if(log.p) -Inf else 0
    out[q >= 1] <- if(log.p) 0 else 1
  } else {
    out[q >= 1] <- dxbeta(1, mu = mu[q >= 1], phi = phi[q >= 1], nu = nu[q >= 1], log = log.p)
    out[q <= 0] <- if(log.p) 0 else 1
    out[q >  1] <- if(log.p) -Inf else 0
  }

  return(out)
}

qxbeta <- function(p, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  q <- qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi, lower.tail = lower.tail, log.p = log.p)
  q <- q * (1 + 2 * nu) - nu
  q[q < 0] <- 0
  q[q > 1] <- 1
  return(q)
}

rxbeta <- function(n, mu, phi, nu = 0) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  r <- -nu + (1 + 2 * nu) * rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  r[r < 0] <- 0
  r[r > 1] <- 1
  return(r)
}

mean_xbeta <- function(mu, phi, nu, ...) {
    a <- mu * phi
    b <- (1 - mu) * phi
    d <- (1 + 2 * nu)
    q0 <- nu / d
    q1 <- (1 + nu) / d
    t3 <- pbeta(q1, a, b)
    t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
    t2 <- nu * (t3 - pbeta(q0, a, b))
    1 + t1 - t2 - t3
}

var_xbeta <- function(mu, phi, nu, quad = 20, ...) {
    if(length(quad) == 1L) quad <- quadtable(quad)
    a <- mu * phi
    b <- (1 - mu) * phi
    mu1 <- (phi * mu + 1) / (phi + 1)
    d <- (1 + 2 * nu)
    q0 <- nu / d
    q1 <- (1 + nu) / d
    t3 <- pbeta(q1, a, b)
    t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
    t2 <- nu * (t3 - pbeta(q0, a, b))
    v1 <- d^2 * mu * mu1 * (pbeta(q1, a + 2, b) - pbeta(q0, a + 2, b))
    v2 <- nu * t2
    v3 <- 2 * nu * t1
    out <- c(v1 + v2 - v3 - t3, t1 - t2 - t3)
    out[1] - out[2] * (2 + out[2])
}


## distributions3 interface

XBeta <- function(mu, phi, nu = 0) {
  n <- c(length(mu), length(phi), length(nu))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  d <- data.frame(mu = mu, phi = phi, nu = nu)
  class(d) <- c("XBeta", "distribution")
  d
}

mean.XBeta <- function(x, ...) {
  m <- vapply(seq_along(x), function(i) mean_xbeta(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(m, names(x))
}

variance.XBeta <- function(x, ...) {
  v <- vapply(seq_along(x), function(i) var_xbeta(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(v, names(x))
}

skewness.XBeta <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.XBeta <- function(x, ...) {
  stop("not yet implemented")
}

random.XBeta <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rxbeta(n = at, mu = d$mu, phi = d$phi, nu = d$nu)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.XBeta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbeta(x = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.XBeta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbeta(x = at, mu = d$mu, phi = d$phi, nu = d$nu, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.XBeta <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pxbeta(q = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.XBeta <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qxbeta(p = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.XBeta <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(rep.int(0, length(d)), rep.int(1, length(d)), d, drop = drop)
}

is_discrete.XBeta <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.XBeta <- function(d, ...) {
  setNames(d$nu <= 0, names(d))
}
