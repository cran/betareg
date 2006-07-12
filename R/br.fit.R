"br.fit" <- function (x, y, link)
{
    x <- as.matrix(x)
    y <- as.matrix(y)
    linkfun <- link$linkfun
    linkinv <- link$linkinv
    mu.eta <- link$mu.eta
    diflink <- link$diflink
    T <- link$T
    ynew = linkfun(y)
    ystar = log(y/(1-y))
    ajuste = lm.fit(x, ynew)
    beta = c(ajuste$coef)
    k = length(beta)
    n = length(y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n - k) * (dlink)^2)
    phi = 1/n * sum(mean * (1 - mean)/sigma2 - 1)
    reg = c(beta, phi)
    loglik <- function(z) {
        z1 = z[1:k]
        z2 = z[k + 1]
        mu = linkinv(x %*% z1)
        sum(lgamma(z2) - lgamma(mu * z2) - lgamma((1 - mu) *
            z2) + (mu * z2 - 1) * log(y) + ((1 - mu) * z2 - 1) *
            log(1 - y))
    }
    loglikt <- function(z) {
        d = length(z) - 1
        z1 = z[1:d]
        z2 = z[d + 1]
        mu = z1
        lgamma(z2) - lgamma(mu * z2) - lgamma((1 - mu) * z2) +
            (mu * z2 - 1) * log(y) + ((1 - mu) * z2 - 1) * log(1 -
            y)
    }
    escore <- function(z) {
        z1 = z[1:k]
        z2 = z[k + 1]
        mu = linkinv(x %*% z1)
        munew = digamma(mu * z2) - digamma((1 - mu) * z2)
        T = diag(c(mu.eta(x %*% z1) ) )
        c(z2 * t(x) %*% T %*% (ystar - munew), sum(digamma(z2) -
            mu * digamma(mu * z2) - (1 - mu) * digamma((1 - mu) *
            z2) + mu * log(y) + (1 - mu) * log(1 - y)))
    }
    opt <- optim(reg, loglik, escore, method = "BFGS", control = list(fnscale = -1,
        maxit = 2000, reltol = 1e-12))
    if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")
    z <- c()
    coef <- (opt$par)[1:k]
    z$coeff <- coef
    z$beta <- beta
    z$phi <- phi
    etahat <- x %*% coef
    phihat <- opt$par[k + 1]
    muhat = linkinv(etahat)
    z$fitted <- muhat
    z$phiest <- phihat
    psi1 = trigamma(muhat * phihat)
    psi2 = trigamma((1 - muhat) * phihat)
    T1 = T(etahat)


    x1 <- rep(1,length(ynew))
    x1 <- as.matrix(x1)
    loglik2 <- function(z) {
        mu = linkinv(x1 %*% z)
        sum(lgamma(phihat) - lgamma(mu * phihat) - lgamma((1 - mu) *
            phihat) + (mu * phihat - 1) * log(y) + ((1 - mu) * phihat - 1) *
            log(1 - y))
    }    
    escore2 <- function(z) {
        mu = linkinv(x1 %*% z)
        munew = digamma(mu * phihat) - digamma((1 - mu) * phihat)
        T = diag(c(mu.eta(x1 %*% z) ) )
        phihat * t(x1) %*% T %*% (ystar - munew)
    }
    ajuste2 = lm.fit(x1, ynew)
    beta1 <- ajuste2$coef
    opt2 <- optim(beta1, loglik2, escore2, method = "BFGS", control = list(fnscale = -1,
        maxit = 2000, reltol = 1e-12))
    if (opt$conv != 0)
        warning("FUNCTION DID NOT CONVERGE!")    
    coef0 <- opt2$par
    etahat0 <- x1%*%coef0
    mu0 <- linkinv(etahat0)
    val2 <- 2 * (loglikt(c(y, phihat)) - loglikt(c(mu0, phihat)))
val2[val2<0] <- 0
    resd2 <- sign(y - mu0) * sqrt(val2)
    nulldev <- sum(resd2^2)
    z$nulldev <- nulldev



    W = diag(c(phihat * (psi1 + psi2))) %*% T1^2
    vc = phihat * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(phihat)))
    tempinv = solve(t(x) %*% W %*% x)
    g = sum(diag(D)) - (1/phihat) * t(vc) %*% t(T1) %*% x %*%
        tempinv %*% t(x) %*% T1 %*% vc
    K1 = tempinv %*% (c(g) * diag(k) + (1/phihat) * t(x) %*%
        T1 %*% vc %*% t(vc) %*% t(T1) %*% x %*% tempinv)
    K2 = -tempinv %*% t(x) %*% T1 %*% vc
    tempmatrix <- (-t(vc) %*% t(T1) %*% x %*% tempinv)
    tempmatrix <- cbind(tempmatrix, phihat)
    fisherinv = (1/(phihat * c(g))) * rbind(cbind(K1, K2), tempmatrix)
    stderrors <- sqrt(diag(fisherinv))[1:k]
    z$stderrors <- stderrors
    phier <- sqrt(diag(fisherinv))[k + 1]
    muhat <- as.vector(muhat)
    H = sqrt(W) %*% x %*% tempinv %*% t(x) %*% sqrt(W)
    h = diag(H)
    z$k = k
    z$h = h
    mustar = digamma(muhat * phihat) - digamma((1 - muhat) *
        phihat)
    Q = (phihat * (trigamma(muhat * phihat) + trigamma((1 - muhat) *
        phihat)) - (ystar - mustar) * etahat/(mu.eta(etahat))) *(mu.eta(etahat))^2
    Q <- as.vector(Q)
    Q <- diag(Q)
    f = vc - (ystar - mustar)
    e = -(y - muhat)/(y * (1 - y))
    XQXinv = solve(t(x) %*% Q %*% x)
    M = 1/(y * (1 - y))
    M = as.vector(M)
    M = diag(M)
    g = sum(diag(D)) - (1/phihat) * t(f) %*% t(T1) %*% x %*% tempinv %*% t(x) %*% T1 %*% f
    GL1 = T1 %*% x %*% XQXinv %*% t(x) %*% T1 %*% M
    GL2 = (1/(c(g) * phihat)) * T1 %*% x %*% XQXinv %*% t(x) %*%
        T1 %*% f %*% (t(f) %*% T1 %*% x %*% XQXinv %*% t(x) %*%
        T1 %*% M - t(e))
    GL = GL1 + GL2
    z$GL = GL
    z$fitted <- muhat
val <- 2 * (loglikt(c(y, phihat)) - loglikt(c(muhat, phihat)))
val[val<0] <- 0
    resd <- sign(y - muhat) * sqrt(val)
resstd <- sqrt(1.0+phihat)*(y-muhat)/sqrt((1.0-h)*muhat*(1-muhat));
z$resstd <- resstd
    z$resd <- resd
    z$phistd <- phier
    z$value <- opt$value

    z$zstats <- coef/stderrors
    z$linpred <- etahat
    res <- y - muhat
    res <- as.vector(res)
    z$res <- res
    z$df.residual <- n-k
    z$pvalues <- 2 * (1 - pnorm(abs(coef/stderrors)))
    if(!(var(etahat)*var(ynew) == 0)){
       pseudor2 <- cor(etahat, ynew)^2       
    }
    else {pseudor2 <- NA}
       z$pseudor2 <- pseudor2
       z$etahat <- etahat
       sigma2 = sum(res^2)/((n - k) * (diflink(muhat)^2))
       z$sigma2 <- sigma2
    z
}