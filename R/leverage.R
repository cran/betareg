## different leverage measures
## both currently have implementation that require
## n x n computations...hence these computations are done
## in the extractor functions (not the fitting function)

gleverage <- function(model, ...) {
  UseMethod("gleverage")
}

gleverage.betareg <- function(model, ...)
{
  ## extract response y and regressors X
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model, model = "mean") else model$x$mean
  z <- if(is.null(model$x)) model.matrix(model, model = "precision") else model$x$precision
  offset <- if(is.null(model$offset)) rep(0, NROW(x)) else model$offset
  wts <- weights(model)
  if(is.null(wts)) wts <- 1
  ystar <- qlogis(y)
  
  ## extract coefficients
  beta <- model$coefficients$mean
  gamma <- model$coefficients$precision

  ## compute different types of "fitted" values
  eta <- as.vector(x %*% beta + offset)
  phi_eta <- as.vector(z %*% gamma)
  mu <- model$link$mean$linkinv(eta)
  phi <- model$link$precision$linkinv(phi_eta)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  #FIXME# Correct generalization??

  ## compute diagonal of T
  Tdiag <- as.vector(model$link$mean$mu.eta(eta))

  ## compute w
  w <- wts * phi * (psi1 + psi2) * Tdiag^2
  ## compute vector c
  vc <- phi * (psi1 * mu - psi2 * (1 - mu))  
  ## compute d
  d <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
  ## compute (X'W X)^(-1)
  xwx1 <- chol2inv(qr.R(qr(sqrt(w) * x)))
  ## compute X'Tc
  xtc <- as.vector(t(x) %*% (Tdiag * vc))  
  ## compute gamma
  xwtc <- as.vector(xwx1 %*% xtc)
  gamma <- sum(d) - sum(xtc * xwtc)/phi  

  ## compute q, f, m, b
  q <- w + wts * (ystar - mustar) * eta/model$link$mean$mu.eta(eta) * Tdiag^2
  ##     -   (in previous version)
  f <- vc - (ystar - mustar)
  m <- 1 / (y * (1 - y))
  b <- -(y - mu) * m

  ## compute T X (X'Q X)^(-1) X' T
  tx_xwx1_xt <- t(t(Tdiag * (x %*% chol2inv(qr.R(qr(sqrt(q) * x))) %*% t(x))) * Tdiag)
  
  ## GL
  gl <- diag(tx_xwx1_xt) * m
  gl <- gl + diag(tx_xwx1_xt %*% f %*% (t(f) %*% tx_xwx1_xt %*% diag(m) - t(b)))/(gamma * phi)
  gl
}

hatvalues.betareg <- function(model, ...)
{
  ## extract response y and regressors X and Z
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model, model = "mean") else model$x$mean
  z <- if(is.null(model$x)) model.matrix(model, model = "precision") else model$x$precision
  offset <- if(is.null(model$offset)) rep(0, NROW(x)) else model$offset
  wts <- weights(model)
  if(is.null(wts)) wts <- 1
  
  ## extract coefficients
  beta <- model$coefficients$mean
  gamma <- model$coefficients$precision

  ## compute different types of "fitted" values
  eta <- as.vector(x %*% beta + offset)
  phi_eta <- as.vector(z %*% gamma)
  mu <- model$link$mean$linkinv(eta)
  phi <- model$link$precision$linkinv(phi_eta)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)

  ## compute w
  w <- wts * phi * (psi1 + psi2) * as.vector(model$link$mean$mu.eta(eta))^2

  ## compute (X'W X)^(-1)
  xwx1 <- chol2inv(qr.R(qr(sqrt(w) * x)))

  ## hat values
  ## (computations are bad: of order n x n)
  w * diag(x %*% xwx1 %*% t(x))
}
