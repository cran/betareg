## four-parameter beta distribution in regression parametrization
## (mean = mu, precision = phi, support = (theta1, theta2))

dbeta4 <- function(x, mu, phi, theta1 = 0, theta2 = 1 - theta1, log = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= theta1 & mu <= theta2),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'theta1' must always be less than 'theta2'" = all(theta1 <= theta2)
  )
  out <- dbeta((x - theta1) / (theta2 - theta1),
    shape1 = mu * phi, shape2 = (1 - mu) * phi, log = log
  )
  out[x <= theta1 | x >= theta2] <- if(log) -Inf else 0
  out <- if(log) out - log(theta2 - theta1) else out/(theta2 - theta1)
  return(out)
}

pbeta4 <- function(q, mu, phi, theta1 = 0, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= theta1 & mu <= theta2),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'theta1' must always be less than 'theta2'" = all(theta1 <= theta2)
  )
  out <- pbeta((q - theta1) / (theta2 - theta1),
    shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p
  )
  return(out)
}

qbeta4 <- function(p, mu, phi, theta1 = 0, theta2 = 1 - theta1, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= theta1 & mu <= theta2),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'theta1' must always be less than 'theta2'" = all(theta1 <= theta2)
  )
  q <- qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi,
    lower.tail = lower.tail, log.p = log.p)
  q <- q * (theta2 - theta1) + theta1
  return(q)
}

rbeta4 <- function(n, mu, phi, theta1 = 0, theta2 = 1 - theta1) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= theta1 & mu <= theta2),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'theta1' must always be less than 'theta2'" = all(theta1 <= theta2)
  )
  r <- theta1 + (theta2 - theta1) * rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  return(r)
}


## distributions3 interface for regression specification
## (mean = mu, precision = phi, support = (theta1, theta2), ...)

Beta4 <- function(mu, phi, theta1 = 0, theta2 = 1 - theta1) {
  n <- c(length(mu), length(phi), length(theta1), length(theta2))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))

  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= theta1 & mu <= theta2),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'theta1' must always be less than 'theta2'" = all(theta1 <= theta2)
  )

  d <- data.frame(mu = mu, phi = phi, theta1 = theta1, theta2 = theta2)
  class(d) <- c("Beta4", "distribution")
  d
}

mean.Beta4 <- function(x, ...) {
  rval <- x$theta1 + (x$theta2 - x$theta1) * x$mu
  setNames(rval, names(x))
}

variance.Beta4 <- function(x, ...) {
  rval <- (x$theta2 - x$theta1)^2 * x$mu * (1 - x$mu)/(1 + x$phi)
  setNames(rval, names(x))
}

skewness.Beta4 <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.Beta4 <- function(x, ...) {
  stop("not yet implemented")
}

random.Beta4 <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rbeta4(n = at, mu = d$mu, phi = d$phi, theta1 = d$theta1, theta2 = d$theta2)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.Beta4 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbeta4(x = at, mu = d$mu, phi = d$phi, theta1 = d$theta1, theta2 = d$theta2, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.Beta4 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbeta4(x = at, mu = d$mu, phi = d$phi, theta1 = d$theta1, theta2 = d$theta2, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.Beta4 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pbeta4(q = at, mu = d$mu, phi = d$phi, theta1 = d$theta1, theta2 = d$theta2, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.Beta4 <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qbeta4(p = at, mu = d$mu, phi = d$phi, theta1 = d$theta1, theta2 = d$theta2, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.Beta4 <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(d$theta1, d$theta2, d, drop = drop)
}

is_discrete.Beta4 <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.Beta4 <- function(d, ...) {
  setNames(rep.int(TRUE, length(d)), names(d))
}
