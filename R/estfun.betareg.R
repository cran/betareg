estfun.betareg <- function(x, ...)
{
  ## extract response y and regressors X
  xmat <- x$x
  y <- x$y
  
  ## extract coefficients
  beta <- coef(x)
  phi <- beta[length(beta)]
  beta <- beta[-length(beta)]

  ## compute y*
  ystar = x$funlink$linkfun(y)
  
  ## compute mu*
  eta <- xmat %*% beta
  mu <- x$linkinv(eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ystar <- as.vector(ystar)
  mustar <- as.vector(mustar)  

  ## compute diagonal of matrix T
  Tdiag <- x$funlink$mu.eta(eta)

  Tdiag <- as.vector(Tdiag)

  ## compute scores of beta
  rval <- phi * (ystar - mustar) * Tdiag * xmat

  ## combine with scores of phi
  rval <- cbind(rval,
    phi = (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)))

  attr(rval, "assign") <- NULL
  return(rval)
}