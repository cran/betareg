"br.fit" <-
function (x, y, link) 
{
    x <- as.matrix(x)
    y <- as.matrix(y)
    linkfun <- link$linkfun
    linkinv <- link$linkinv
    mu.eta <- link$mu.eta
    diflink <- link$diflink
    T <- link$T
    ynew = linkfun(y)
    ajuste = lm.fit(x, ynew)
    beta = c(ajuste$coef)
    k = length(beta)
    n = length(y)
    mean = fitted(ajuste)
    mean = exp(mean)/(1 + exp(mean))
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n - k) * (dlink)^2)
    phi = 1/n * sum(mean * (1 - mean)/sigma2 - 1)
    reg = c(beta, phi)
    loglik <- function(z) {
        z1 = z[1:k]
        z2 = z[k+1]
        mu = linkinv(x%*%z1)
        sum(lgamma(z2) - lgamma(mu * z2) - lgamma((1 - 
            mu) * z2) + (mu * z2 - 1) * log(y) + ((1 - mu) * 
            z2 - 1) * log(1 - y))
    }
    loglikt <- function(z){
        d = length(z)-1
        z1 = z[1:d]
        z2 = z[d+1]
        mu = z1
        lgamma(z2) - lgamma(mu * z2) - lgamma((1 - 
            mu) * z2) + (mu * z2 - 1) * log(y) + ((1 - mu) * 
            z2 - 1) * log(1 - y)
    }
    escore <- function(z) {
        z1 = z[1:k]
        z2 = z[k + 1]
        mu = linkinv(x %*% z1)
        munew = digamma(mu * z2) - digamma((1 - mu) * z2)
        T = diag(c(exp(x %*% z1)/(1 + exp(x %*% z1))^2))
        c(z2 * t(x) %*% T %*% (ynew - munew), sum(digamma(z2) - 
            mu * digamma(mu * z2) - (1 - mu) * digamma((1 - mu) * 
            z2) + mu * log(y) + (1 - mu) * log(1 - y)))
    }
    opt <- optim(reg, loglik, escore, method = "BFGS", control = list(fnscale = -1,maxit=2000))
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
    H = sqrt(W)%*%x%*%tempinv%*%t(x)%*%sqrt(W)
    h = diag(H)
    z$k = k
    z$h = h
    ystar = ynew 
    mustar = digamma(muhat*phihat) - digamma((1.0-muhat)*phihat)
    Q = (phihat*(trigamma(muhat*phihat) + trigamma((1-muhat)*phihat)) - (ystar-mustar)*(1-2*muhat)/(muhat*(1-muhat)))*(muhat^2)*(1-muhat)^2
    Q <- as.vector(Q)
    Q <- diag(Q)
    f = vc - (ystar-mustar)
    e = -(y-muhat)/(y*(1-y))
    XQXinv = solve(t(x)%*%Q%*%x)
    M = 1/(y*(1-y))
    M = as.vector(M)
    M = diag(M)
    GL1 = T1%*%x%*%XQXinv%*%t(x)%*%T1%*%M
    GL2 = (1/(c(g)*phihat))*T1%*%x%*%XQXinv%*%t(x)%*%T1%*%f%*%(t(f)%*%T1%*%x%*%XQXinv%*%t(x)%*%T1%*%M-t(e))
    GL = GL1 + GL2
    z$GL = GL
    z$fitted <- muhat
    resd<-sign(y-muhat)*sqrt(2*(loglikt(c(y,phihat)) - loglikt(c(muhat,phihat))))
    z$resd <- resd
    z$phistd <- phier
    z$zstats <- coef/stderrors
    res <- y - muhat
    res <- as.vector(res)
    z$res <- res
    z$pvalues <- 2 * (1 - pnorm(abs(coef/stderrors)))
    pseudor2 <- cor(etahat, ynew)^2
    z$pseudor2 <- pseudor2
    z
}
