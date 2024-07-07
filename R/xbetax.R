## extended-domain beta mixture distribution (XBetaX)
## based exponential mixture of the censored symmetric four-parameter beta distribution in regression parameterization
## (mean = mu, precision = phi, latent support = (-Inf, Inf) censored to [0, 1])

## auxiliary quadrature function
quadtable <- function(nquad = 20) {
  matrix(unlist(
    statmod::gauss.quad(n = nquad, kind = "laguerre", alpha = 0)
  ), nrow = nquad)
}

dxbetax <- function(x, mu, phi, nu = 0, log = FALSE, quad = 20) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) return(dbeta(x, shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log))

  ## unify lengths of all variables
  n <- max(length(x), length(mu), length(phi), length(nu))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  if(length(quad) == 1L) quad <- quadtable(quad)
  out <- apply(quad, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * dbeta((x + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)/(1 + 2 * e)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  ## censoring
  out[x <= 0] <- pxbetax(0, mu = mu[x <= 0], phi = phi[x <= 0], nu = nu[x <= 0])
  out[x >= 1] <- pxbetax(1, mu = mu[x >= 1], phi = phi[x >= 1], nu = nu[x >= 1], lower.tail = FALSE)
  out[x < 0 | x > 1] <- 0

  ## additional arguments
  if(log) out <- log(out)

  return(out)
}

pxbetax <- function(q, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) return(pbeta(q, shape1 = mu * phi, shape2 = (1 - mu) * phi, lower.tail = lower.tail, log.p = log.p))

  ## unify lengths of all variables
  n <- max(length(q), length(mu), length(phi), length(nu))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)

  ## quadrature
  if(length(quad) == 1L) quad <- quadtable(quad)
  out <- apply(quad, 1, function(rule) {
    e <- rule[1] * nu
    rule[2] * pbeta((q + e)/(1 + 2 * e), shape1 = mu * phi, shape2 = (1 - mu) * phi)
  })
  out <- if (is.null(dim(out))) sum(out) else rowSums(out)

  ## tail of the distribution
  if(!lower.tail) out <- 1 - out

  ## censoring
  if(lower.tail) {
    out[q <  0] <- 0
    out[q >= 1] <- 1
    out[q <= 0 & nu <= 0] <- 0 ## for nu = 0 no point mass at 0
  } else {
    out[q <= 0] <- 1
    out[q >  1] <- 0
    out[q >= 1 & nu <= 0] <- 0 ## for nu = 0 no point mass at 1
  }

  ## additional arguments
  if(log.p) out <- log(out)
  return(out)
}

qxbetax <- function(p, mu, phi, nu = 0, lower.tail = TRUE, log.p = FALSE, quad = 20, tol = .Machine$double.eps^0.7) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )

  ## unify lengths of all variables
  n <- max(length(p), length(mu), length(phi), length(nu))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  nu <- rep_len(nu, n)
  q <- rep_len(NA_real_, n)

  ## quadrature
  if(length(quad) == 1L) quad <- quadtable(quad)

  ## cumulative probabilities at boundary
  p0 <- pxbetax(0, mu = mu, phi = phi, nu = nu, log.p = log.p, quad = quad, lower.tail = TRUE)
  p1 <- pxbetax(1, mu = mu, phi = phi, nu = nu, log.p = log.p, quad = quad, lower.tail = FALSE)
  p1 <- if(!log.p) 1 - p1 else log1p(-exp(p1))

  ## indexes for boundary and non-boundary observations
  idx0 <- if (lower.tail) p <= p0 else if (!log.p) 1 - p <= p0 else log1p(-exp(p)) <= p0
  idx1 <- if (lower.tail) p >= p1 else if (!log.p) 1 - p >= p1 else log1p(-exp(p)) >= p1
  idx <- !idx0 & !idx1

  ## boundary quantiles
  if (any(idx0)) q[idx0] <- 0
  if (any(idx1)) q[idx1] <- 1

  ## non-boundary quantiles
  if (any(idx)) {
    obj <- function(pq, mu, phi, nu, p) {
      p - pxbetax(q = pq, mu = mu, phi = phi, nu = nu, lower.tail = lower.tail, log.p = log.p, quad = quad)
    }
    iroot <- function(i) {
      r <- try(uniroot(obj, c(0, 1), mu = mu[i], phi = phi[i], nu = nu[i], p = p[i], tol = tol), silent = TRUE)
      if(inherits(r, "try-error")) NA_real_ else r$root
    }
    q[idx] <- vapply(which(idx), iroot, 0.0)
  }

  return(q)
}

rxbetax <- function(n, mu, phi, nu = 0) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  if(isTRUE(all(nu == 0))) {
    rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  } else {
    rxbeta(n, mu = mu, phi = phi, nu = rexp(n, 1/nu))
  }
}


mean_xbetax <- function(mu, phi, nu, quad = 20, ...) {
    if(length(quad) == 1L) quad <- quadtable(quad)
    a <- mu * phi
    b <- (1 - mu) * phi
    out <- apply(quad, 1, function(rule) {
        e <- rule[1] * nu
        d <- (1 + 2 * e)
        q0 <- e / d
        q1 <- (1 + e) / d
        t3 <- pbeta(q1, a, b)
        t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
        t2 <- e * (t3 - pbeta(q0, a, b))
        rule[2] * (t1 - t2 - t3)
    })
    1 + sum(out)
}

var_xbetax <- function(mu, phi, nu, quad = 20, ...) {
    if(length(quad) == 1L) quad <- quadtable(quad)
    a <- mu * phi
    b <- (1 - mu) * phi
    mu1 <- (phi * mu + 1) / (phi + 1)
    out <- apply(quad, 1, function(rule) {
        e <- rule[1] * nu
        d <- (1 + 2 * e)
        q0 <- e / d
        q1 <- (1 + e) / d
        t3 <- pbeta(q1, a, b)
        t1 <- d * mu * (pbeta(q1, a + 1, b) - pbeta(q0, a + 1, b))
        t2 <- e * (t3 - pbeta(q0, a, b))
        v1 <- d^2 * mu * mu1 * (pbeta(q1, a + 2, b) - pbeta(q0, a + 2, b))
        v2 <- e * t2
        v3 <- 2 * e * t1
        rule[2] * c(v1 + v2 - v3 - t3, t1 - t2 - t3)
    })
    out <- rowSums(out)
    out[1] - out[2] * (2 + out[2])
}



## distributions3 interface

XBetaX <- function(mu, phi, nu = 0) {
  n <- c(length(mu), length(phi), length(nu))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'nu' must always be non-negative" = all(nu >= 0)
  )
  d <- data.frame(mu = mu, phi = phi, nu = nu)
  class(d) <- c("XBetaX", "distribution")
  d
}

mean.XBetaX <- function(x, ...) {
  m <- vapply(seq_along(x), function(i) mean_xbetax(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(m, names(x))
}

variance.XBetaX <- function(x, ...) {
  v <- vapply(seq_along(x), function(i) var_xbetax(mu = x$mu[i], phi = x$phi[i], nu = x$nu[i], ...), 0.0)
  setNames(v, names(x))
}

skewness.XBetaX <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.XBetaX <- function(x, ...) {
  stop("not yet implemented")
}

random.XBetaX <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rxbetax(n = at, mu = d$mu, phi = d$phi, nu = d$nu)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbetax(x = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dxbetax(x = at, mu = d$mu, phi = d$phi, nu = d$nu, log = TRUE, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.XBetaX <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pxbetax(q = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.XBetaX <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qxbetax(p = at, mu = d$mu, phi = d$phi, nu = d$nu, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.XBetaX <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(rep.int(0, length(d)), rep.int(1, length(d)), d, drop = drop)
}

is_discrete.XBetaX <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.XBetaX <- function(d, ...) {
  setNames(d$nu <= 0, names(d))
}
