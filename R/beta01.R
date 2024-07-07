## zero- and/or one-inflated beta distribution in regression parametrization
## (mean = mu, precision = phi, p0 = probability for 0, p1 = probability for 1)

dbeta01 <- function(x, mu, phi, p0 = 0, p1 = 0, log = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'p0' must always be in [0, 1]" = all(p0 >= 0 & p0 <= 1),
    "parameter 'p1' must always be in [0, 1]" = all(p1 >= 0 & p1 <= 1),
    "sum of parameters 'p0' + 'p1' must always be in [0, 1]" = all(p0 + p1 <= 1)
  )
  rval <- dbetar(x, mu = mu, phi = phi, log = log)

  ## unify lengths of variables
  n <- length(rval)
  x <- rep_len(x, n)
  p0 <- rep_len(p0, n)
  p1 <- rep_len(p1, n)
  
  if(log) {
    rval <- rval + log(1 - p0 - p1)
    rval[x <= 0] <- log(p0[x <= 0])
    rval[x >= 1] <- log(p1[x >= 1])
    rval[x < 0 | x > 1] <- -Inf
  } else {
    rval <- rval * (1 - p0 - p1)
    rval[x <= 0] <- p0[x <= 0]
    rval[x >= 1] <- p1[x >= 1]
    rval[x < 0 | x > 1] <- 0
  }
  return(rval)
}

pbeta01 <- function(q, mu, phi, p0 = 0, p1 = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'p0' must always be in [0, 1]" = all(p0 >= 0 & p0 <= 1),
    "parameter 'p1' must always be in [0, 1]" = all(p1 >= 0 & p1 <= 1),
    "sum of parameters 'p0' + 'p1' must always be in [0, 1]" = all(p0 + p1 <= 1)
  )
  rval <- p0 + (1 - p0 - p1) * pbetar(q, mu = mu, phi = phi)

  ## unify lengths of variables
  n <- length(rval)
  q <- rep_len(q, n)
  p0 <- rep_len(p0, n)
  p1 <- rep_len(p1, n)

  if(lower.tail) { 
    rval[q <= 0] <- p0[q <= 0]
    rval[q <  0] <- 0
    rval[q >= 1] <- 1
  } else {
    rval <- 1 - rval
    rval[q >= 1] <- p1[q >= 1]
    rval[q >  1] <- 0
    rval[q <= 0] <- 1
  }
  if(log.p) rval <- log(rval)
  return(rval)
}

qbeta01 <- function(p, mu, phi, p0 = 0, p1 = 0, lower.tail = TRUE, log.p = FALSE) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'p0' must always be in [0, 1]" = all(p0 >= 0 & p0 <= 1),
    "parameter 'p1' must always be in [0, 1]" = all(p1 >= 0 & p1 <= 1),
    "sum of parameters 'p0' + 'p1' must always be in [0, 1]" = all(p0 + p1 <= 1)
  )
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  rval <- qbetar(pmin(pmax((p - p0)/(1 - p0 - p1), 0), 1), mu = mu, phi = phi)
  n <- length(rval)
  p <- rep_len(p, n)
  rval[p <= p0] <- 0
  rval[p > (1 - p1)] <- 1
  rval[p0 < 0 | p0 > 1 | p1 < 0 | p1 > 1] <- NaN
  return(rval)
}

rbeta01 <- function(n, mu, phi, p0 = 0, p1 = 0) {
  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'p0' must always be in [0, 1]" = all(p0 >= 0 & p0 <= 1),
    "parameter 'p1' must always be in [0, 1]" = all(p1 >= 0 & p1 <= 1),
    "sum of parameters 'p0' + 'p1' must always be in [0, 1]" = all(p0 + p1 <= 1)
  )

  ## unify lengths of variables
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  p0 <- rep_len(p0, n)
  p1 <- rep_len(p1, n)

  ## sample left boundary, right boundary, middle
  rval <- rbinom(n, size = 1, prob = 1 - p0)
  i <- rval > 0
  rval[i] <- rval[i] + rbinom(sum(i), size = 1, prob = (1 - p0 - p1)/(1 - p0))
  i <- rval > 1
  rval[i] <- rbetar(sum(i), mu = mu[i], phi = phi[i])
  return(rval)
}


## distributions3 interface

Beta01 <- function(mu, phi, p0 = 0, p1 = 0) {
  n <- c(length(mu), length(phi), length(p0), length(p1))
  stopifnot("parameter lengths do not match (only scalars are allowed to be recycled)" = all(n %in% c(1L, max(n))))

  stopifnot(
    "parameter 'mu' must always be in [0, 1]" = all(mu >= 0 & mu <= 1),
    "parameter 'phi' must always be non-negative" = all(phi >= 0),
    "parameter 'p0' must always be in [0, 1]" = all(p0 >= 0 & p0 <= 1),
    "parameter 'p1' must always be in [0, 1]" = all(p1 >= 0 & p1 <= 1),
    "sum of parameters 'p0' + 'p1' must always be in [0, 1]" = all(p0 + p1 <= 1)
  )

  d <- data.frame(mu = mu, phi = phi, p0 = p0, p1 = p1)
  class(d) <- c("Beta01", "distribution")
  d
}

mean.Beta01 <- function(x, ...) {
  rval <- (1 - x$p0 - x$p1) * x$mu + x$p1
  setNames(rval, names(x))
}

variance.Beta01 <- function(x, ...) {
  pm <- 1 - x$p0 - x$p1
  rval <- x$p1 * (1 - x$p1) + x$mu^2 * pm * (1 - pm) - 2 * x$mu * pm * x$p1 + pm * x$mu * (1 - x$mu)/(1 + x$phi)
  setNames(rval, names(x))
}

skewness.Beta01 <- function(x, ...) {
  stop("not yet implemented")
}

kurtosis.Beta01 <- function(x, ...) {
  stop("not yet implemented")
}

random.Beta01 <- function(x, n = 1L, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  n <- distributions3::make_positive_integer(n)
  if (n == 0L) return(numeric(0L))
  FUN <- function(at, d) rbeta01(n = at, mu = d$mu, phi = d$phi, p0 = d$p0, p1 = d$p1)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = n, type = "random", drop = drop)
}

pdf.Beta01 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbeta01(x = at, mu = d$mu, phi = d$phi, p0 = d$p0, p1 = d$p1, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "density", drop = drop, elementwise = elementwise)
}

log_pdf.Beta01 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) dbeta01(x = at, mu = d$mu, phi = d$phi, p0 = d$p0, p1 = d$p1, log = TRUE)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "logLik", drop = drop, elementwise = elementwise)
}

cdf.Beta01 <- function(d, x, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) pbeta01(q = at, mu = d$mu, phi = d$phi, p0 = d$p0, p1 = d$p1, ...)
  distributions3::apply_dpqr(d = d, FUN = FUN, at = x, type = "probability", drop = drop, elementwise = elementwise)
}

quantile.Beta01 <- function(x, probs, drop = TRUE, elementwise = NULL, ...) {
  stopifnot(requireNamespace("distributions3"))
  FUN <- function(at, d) qbeta01(p = at, mu = d$mu, phi = d$phi, p0 = d$p0, p1 = d$p1, ...)
  distributions3::apply_dpqr(d = x, FUN = FUN, at = probs, type = "quantile", drop = drop, elementwise = elementwise)
}

support.Beta01 <- function(d, drop = TRUE, ...) {
  stopifnot(requireNamespace("distributions3"))
  distributions3::make_support(rep.int(0, length(d)), rep.int(1, length(d)), d, drop = drop)
}

is_discrete.Beta01 <- function(d, ...) {
  setNames(rep.int(FALSE, length(d)), names(d))
}

is_continuous.Beta01 <- function(d, ...) {
  setNames(!(d$p0 > 0 | d$p1 > 0), names(d))
}
