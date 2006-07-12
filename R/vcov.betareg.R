vcov.betareg <- function(object, ...)
{
  ## eobjecttract response y and regressors X
  xmat <- object$x
  y <- object$y
  
  ## extract coefficients
  beta <- coef(object)
  k <- length(beta) - 1
  phi <- beta[k+1]
  beta <- beta[-(k+1)]
  
  ## auxiliary quantities
  eta <- xmat %*% beta
  mu <- object$linkinv(eta)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)

  ## compute diagonal of T
  Tdiag <- object$funlink$mu.eta(eta)

  ## compute w
  w <- phi * (psi1 + psi2) * Tdiag^2
  
  ## compute vector c
  vc <- phi * (psi1 * mu - psi2 * (1 - mu))
  
  ## compute d
  d <- psi1*mu^2 + psi2*(1-mu)^2 - trigamma(phi)
  w <- as.vector(w)
  T <- diag(Tdiag)
  W <- diag(w)

  ## compute (X'W X)^(-1)
  xwx1 <- chol2inv(qr.R(qr(sqrt(w)* xmat)))
  #xwx1 <- solve(t(xmat)%*%W%*%xmat)

  ## compute X'Tc
  xtc <- as.vector(t(xmat) %*% (Tdiag * vc))
  
  ## compute gamma
  xwtc <- as.vector(xwx1 %*% xtc)
  gamma <- sum(d) - sum(xtc * xwtc)/phi
  
  ## compute components of K^(-1)
  Kbb <- (xwx1/phi) %*% (diag(k) + outer(xtc, xwtc)/(gamma*phi))
  Kpp <- (1/gamma)
  Kbp <- -as.vector(xwx1 %*% xtc)/(gamma * phi)

  rval <- rbind(cbind(Kbb, Kbp), c(Kbp, Kpp))
  rownames(rval) <- colnames(rval) <- c(colnames(xmat), "phi")
  return(rval)
}