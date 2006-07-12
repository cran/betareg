hatvalues.betareg <- function(model, ...)
{
    if (!inherits(model, "betareg")) 
    stop("Use only with 'betareg' objects")
    xmat <- model$x
    y <- model$y
  
    ## extract coefficients
    beta <- coef(model)
    k <- length(beta) - 1
    phi <- beta[k+1]
    beta <- beta[-(k+1)]
  
    ## auxiliary quantities
    eta <- xmat %*% beta
    mu <- model$linkinv(eta)
    psi1 <- trigamma(mu * phi)
    psi2 <- trigamma((1 - mu) * phi)

    ## compute diagonal of T
    Tdiag <- model$funlink$mu.eta(eta)

    ## compute w
    w <- phi * (psi1 + psi2) * Tdiag^2
    w <- as.vector(w)
    w <- diag(w)
    hat <- sqrt(w)%*%xmat%*%solve(t(xmat)%*%w%*%xmat)%*%t(xmat)%*%sqrt(w)
    hat <- diag(hat)
    hat
}