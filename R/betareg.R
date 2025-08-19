betareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    link.phi = NULL, type = c("ML", "BC", "BR"),
		    dist = NULL, nu = NULL,
                    control = betareg.control(...),
                    model = TRUE, y = TRUE, x = FALSE,
                    ...)
{
    ## call
    cl <- match.call()
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE

    ## formula
    oformula <- as.formula(formula)
    formula <- as.Formula(formula)
    if(length(formula)[2L] < 2L) {
        formula <- as.Formula(formula(formula), ~ 1)
        simple_formula <- TRUE
    } else {
        if(length(formula)[2L] > 2L) {
            formula <- Formula(formula(formula, rhs = 1L:2L))
            warning("formula must not have more than two RHS parts")
        }
        simple_formula <- FALSE
    }
    mf$formula <- formula

    ## evaluate model.frame
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    ## extract terms, model matrix, response
    mt <- terms(formula, data = data)
    mtX <- terms(formula, data = data, rhs = 1L)
    mtZ <- delete.response(terms(formula, data = data, rhs = 2L))
    Y <- model.response(mf, "numeric")
    X <- model.matrix(mtX, mf)
    Z <- model.matrix(mtZ, mf)

    ## obtain correct subset of predvars/dataClasses to terms
    .add_predvars_and_dataClasses <- function(terms, model.frame) {
      ## original terms
      rval <- terms
      ## terms from model.frame
      nval <- if(inherits(model.frame, "terms")) model.frame else terms(model.frame, dot = control$dot)

      ## associated variable labels
      ovar <- sapply(as.list(attr(rval, "variables")), deparse)[-1]
      nvar <- sapply(as.list(attr(nval, "variables")), deparse)[-1]
      if(!all(ovar %in% nvar)) stop(
        paste("The following terms variables are not part of the model.frame:",
        paste(ovar[!(ovar %in% nvar)], collapse = ", ")))
      ix <- match(ovar, nvar)

      ## subset predvars
      if(!is.null(attr(rval, "predvars")))
        warning("terms already had 'predvars' attribute, now replaced")
      attr(rval, "predvars") <- attr(nval, "predvars")[1L + c(0L, ix)]

      ## subset dataClasses
      if(!is.null(attr(rval, "dataClasses")))
        warning("terms already had 'dataClasses' attribute, now replaced")
      attr(rval, "dataClasses") <- attr(nval, "dataClasses")[ix]

      return(rval)
    }
    mt  <- .add_predvars_and_dataClasses(mt,  mf)
    mtX <- .add_predvars_and_dataClasses(mtX, mf)
    mtZ <- .add_predvars_and_dataClasses(mtZ, mf)

    ## sanity checks
    if(length(Y) < 1) stop("empty model")
    if(any(Y < 0 | Y > 1)) stop("invalid dependent variable, all observations must be in [0, 1]")
    n <- length(Y)

    ## weights
    weights <- model.weights(mf)
    if(is.null(weights)) weights <- 1
    if(length(weights) == 1) weights <- rep.int(weights, n)
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)

    ## offsets
    expand_offset <- function(offset) {
        if(is.null(offset)) offset <- 0
        if(length(offset) == 1) offset <- rep.int(offset, n)
        as.vector(offset)
    }
    ## in mean part of formula
    offsetX <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 1L, terms = TRUE)))
    ## in precision part of formula
    offsetZ <- expand_offset(model.offset(model.part(formula, data = mf, rhs = 2L, terms = TRUE)))
    ## in offset argument (used for mean)
    if(!is.null(cl$offset)) offsetX <- offsetX + expand_offset(mf[, "(offset)"])
    ## collect
    offset <- list(mu = offsetX, phi = offsetZ)

    ## type of estimator and distribution
    type <- match.arg(type, c("ML", "BC", "BR"))
    if(!is.null(dist)) {
      dist <- tolower(as.character(dist))
      if(dist == "xb") dist <- "xbeta"
      if(dist == "xbx") dist <- "xbetax"
      dist <- match.arg(dist, c("beta", "xbeta", "xbetax"))
      if(dist == "beta") {
        if(any((Y <= 0 | Y >= 1)[weights > 0])) stop("dependent variable not suitable for 'beta' distribution, all observations must be in (0, 1)")
      }
      if(dist == "xbeta" && is.null(nu)) {
        warning(sprintf("estimation of 'nu' with 'xbeta' distribution is not feasible, using '%s' instead",
          dist <- if(any((Y <= 0 | Y >= 1)[weights > 0])) "xbetax" else "beta"))
      }
    } else {
      if(is.null(nu)) {
        dist <- if(all((Y > 0 & Y < 1)[weights > 0])) "beta" else "xbetax"
      } else {
        dist <- "xbeta"
      }
    }
    if(dist != "beta" && type != "ML") {
      warning(sprintf("only 'ML' estimation is available for '%s' distribution", dist))
      type <- "ML"
    }

    ## links
    if(is.character(link)) link <- match.arg(link)
    if(is.null(link.phi)) link.phi <- if(simple_formula) "identity" else "log"
    if(is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))

    ## call the actual workhorse: betareg.fit()
    rval <- betareg.fit(X, Y, Z, weights, offset, link, link.phi, type, control, dist, nu)

    ## further model information
    rval$call <- cl
    rval$formula <- oformula
    rval$terms <- list(mu = mtX, phi = mtZ, full = mt)
    rval$levels <- list(mu = .getXlevels(mtX, mf), phi = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
    rval$contrasts <- list(mu = attr(X, "contrasts"), phi = attr(Z, "contrasts"))
    if(model) rval$model <- mf
    if(y) rval$y <- Y
    if(x) rval$x <- list(mu = X, phi = Z)
    if(is.null(rval$dist) || (rval$dist == "beta")) {
        for(n in intersect(names(rval), fix_names_mu_phi)) names(rval[[n]])[1L:2L] <- c("mean", "precision")
    }
    class(rval) <- "betareg"
    return(rval)
}

## internal tools for canonicalizing mean vs. mu and precision vs. phi
fix_names_mu_phi <- c("coefficients", "offset", "link", "terms", "levels", "contrasts", "x")
fix_model_mu_phi <- function(model) {
    model <- sapply(model, function(m) {
        if(m %in% c("mu", "phi", "nu")) return(m)
        match.arg(m, c("mean", "precision", "full"))
    })
    model[model == "mean"] <- "mu"
    model[model == "precision"] <- "phi"
    return(as.vector(model))
}

betareg.control <- function(phi = TRUE, method = "BFGS", maxit = 5000, gradient = NULL, hessian = FALSE, trace = FALSE, start = NULL,
                            fsmaxit = 200, fstol = 1e-8, quad = 20,...)
{
    method <- match.arg(method, c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "nlminb"))
    rval <- list(phi = phi, method = method, maxit = maxit, gradient = gradient, hessian = hessian,
                 trace = trace, start = start, fsmaxit = fsmaxit, fstol = fstol, quad = quad)
    rval <- c(rval, list(...))
    if (method != "nlminb") {
        if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
        rval$fnscale <- -1
        if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
    }
    return(rval)
}

betareg.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
                        link = "logit", link.phi = "log", type = "ML", control = betareg.control(),
                        dist = NULL, nu = NULL)
{
    ## estimation type and distribution:
    if(is.null(dist)) {
      if(is.null(nu)) {
        dist <- if(all(y > 0 & y < 1)) "beta" else "xbetax"
      } else {
        dist <- "xbeta"
      }
    }

    ## only plain ML supported for extended-support distributions
    if(dist != "beta") {
        if(type != "ML") stop(sprintf("only 'ML' estimation implemented for '%s'", dist))
        control$hessian <- TRUE
        control$fsmaxit <- 0
    }
    estnu <- if(dist != "beta" && is.null(nu)) TRUE else FALSE

    ## response and regressor matrix
    n <- NROW(x)
    k <- NCOL(x)
    if(is.null(weights)) weights <- rep.int(1, n)
    nobs <- sum(weights > 0)
    if(is.null(offset)) offset <- rep.int(0, n)
    if(!is.list(offset)) offset <- list(mu = offset, phi = rep.int(0, n))
    if(is.null(z)) {
        m <- 1L
        z <- matrix(1, ncol = m, nrow = n)
        colnames(z) <- "(Intercept)"
        rownames(z) <- rownames(x)
        phi_const <- TRUE
    } else {
        m <- NCOL(z)
        if(m < 1L) stop("dispersion regression needs to have at least one parameter")
        phi_const <- (m == 1L) && isTRUE(all.equal(as.vector(z[, 1L]), rep.int(1, n)))
    }

    ## link processing
    if(is.character(link)) {
        linkstr <- link
        if(linkstr != "loglog") {
            linkobj <- make.link(linkstr)
            ## add d2mu.deta potentially needed for BC/BR
            linkobj$d2mu.deta <- make.d2mu.deta(linkstr)
        } else {
            linkobj <- structure(list(
                linkfun = function(mu) -log(-log(mu)),
                linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
                mu.eta = function(eta) {
                    eta <- pmin(eta, 700)
                    pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
                },
                d2mu.deta = function(eta) pmax(exp(-exp(-eta) - eta) * expm1(-eta), .Machine$double.eps),
                valideta = function(eta) TRUE,
                name = "loglog"
            ), class = "link-glm")
        }
    } else {
        linkobj <- link
        linkstr <- link$name
        if(type != "ML") {
            if(is.null(linkobj$dmu.deta) & is.null(linkobj$d2mu.deta)) {
                warning("link needs to provide d2mu.deta component for BC/BR")
            } else {
                if(is.null(linkobj$d2mu.deta)) linkobj$d2mu.deta <- linkobj$dmu.deta
            }
        }
    }
    linkfun <- linkobj$linkfun
    linkinv <- linkobj$linkinv
    mu.eta <- linkobj$mu.eta
    d2mu.deta <- linkobj$d2mu.deta
    if(is.character(link.phi)) {
        phi_linkstr <- link.phi
        phi_linkobj <- make.link(phi_linkstr)
        phi_linkobj$d2mu.deta <- make.d2mu.deta(phi_linkstr)
    } else {
        phi_linkobj <- link.phi
        phi_linkstr <- link.phi$name
        if(type != "ML") {
            if(is.null(phi_linkobj$dmu.deta) & is.null(phi_linkobj$d2mu.deta)) {
                warning("link.phi needs to provide d2mu.deta component for BC/BR")
            } else {
                if(is.null(phi_linkobj$d2mu.deta)) phi_linkobj$d2mu.deta <- phi_linkobj$dmu.deta
            }
        }
    }
    phi_linkfun <- phi_linkobj$linkfun
    phi_linkinv <- phi_linkobj$linkinv
    phi_mu.eta <- phi_linkobj$mu.eta
    phi_d2mu.deta <- phi_linkobj$d2mu.deta
    ## y* transformation
    ystar <- qlogis(y)

    ## control parameters
    ocontrol <- control
    phi_full <- control$phi
    method <- control$method
    gradient <- control$gradient
    hessian <- control$hessian
    start <- control$start
    fsmaxit <- control$fsmaxit
    fstol <- control$fstol
    quad <- control$quad
    control$phi <- control$method <- control$gradient <- control$hessian <- control$start <- control$fsmaxit <- control$fstol <- control$quad <- NULL
    if(is.null(gradient)) gradient <- (dist == "beta") ## currently use analytical gradients only for classic beta

    ## starting values
    if(is.null(start)) {
        auxreg <- lm.wfit(x, if(dist == "beta") linkfun(y) else linkfun((y * (n - 1) + 0.5)/n), weights, offset = offset[[1L]])
        beta <- auxreg$coefficients
        yhat <- linkinv(auxreg$fitted.values)
        dlink <- 1/mu.eta(linkfun(yhat))
        res <- auxreg$residuals
        res[weights <= 0] <- 0
        sigma2 <- sum(weights * res^2)/((sum(weights) - k) * (dlink)^2)
        phi_y <- weights * yhat * (1 - yhat)/(sum(weights) * sigma2) - 1/n
        phi <- rep.int(0, ncol(z))
        phi[1L] <- suppressWarnings(phi_linkfun(sum(phi_y)))
        ## i.e., start out from the fixed dispersion model as described
        ## in Ferrari & Cribari-Neto (2004) (and differing from Simas et al. 2009)
        ## An alternative would be
        ##   phi <- lm.wfit(z, phi_linkfun(phi_y), weights)$coefficients
        ## but that only works in general if all(phi_y > 0) which is not necessarily
        ## the case.
        ##
        ## Additionally, sum(phi_y) might not even be > 0 which should be caught.
        if(!isTRUE(phi_linkinv(phi[1L]) > 0)) {
            warning("no valid starting value for precision parameter found, using 1 instead")
            phi[1L] <- 1
        }
        start <- list(mu = beta, phi = phi)
        if(estnu) start$nu <- log(mean(y <= 0 | y >= 1))
    }
    if(is.list(start)) start <- do.call("c", start)

    indices01 <- (y > 0) & (y < 1)
    indices0 <- (y <= 0)
    indices1 <- (y >= 1)

    ## various fitted quantities (parameters, linear predictors, etc.)
    fitfun <- function(par, deriv = 0L) {
        beta <- par[seq.int(length.out = k)]
        gamma <- par[seq.int(length.out = m) + k]
        nu <- if(estnu) exp(par[k + m + 1]) else nu
        eta <- as.vector(x %*% beta + offset[[1L]])
        phi_eta <- as.vector(z %*% gamma + offset[[2L]])
        mu <- linkinv(eta)
        phi <- phi_linkinv(phi_eta)
        shape1 <- mu * phi
        shape2 <- (1 - mu) * phi
        if (deriv >= 1L) {
            d1 <- digamma(shape1)
            d2 <- digamma(shape2)
            mustar <- d1 - d2
        }
        else {
            d1 <- d2 <- mustar <- NULL
        }
        if (deriv >= 2L) {
            psi1 <- trigamma(shape1)
            psi2 <- trigamma(shape2)
        }
        else {
            psi1 <- psi2 <- NULL
        }
        if (deriv >= 3L) {
            b12 <- beta(shape1, shape2)
            d12 <- digamma(phi)
        }
        else {
            b12 <- d12 <- NULL
        }
        list(
            shape1 = shape1,
            shape2 = shape2,
            d1 = d1,
            d2 = d2,
            b12 = b12,
            d12 = d12,
            beta = beta,
            gamma = gamma,
            nu = nu,
            eta = eta,
            phi_eta = phi_eta,
            mu = mu,
            phi = phi,
            mustar = mustar,
            psi1 = psi1,
            psi2 = psi2
        )
    }

    if(dist == "beta") {

        ## objective function
        loglikfun <- function(par, fit = NULL) {
            ## extract fitted parameters
            if(is.null(fit)) fit <- fitfun(par)
            ## compute log-likelihood
            ## conditionals are still here for backward compatibility; to be
            ## removed in a later release
            with(fit, {
                if(any(!is.finite(phi)) | any(shape1 > 1e300) | any(shape2 > 1e300)) {
                    NaN
                }
                ## catch extreme cases
                else {
                    ll <- suppressWarnings(dbeta(y, shape1, shape2, log = TRUE))
                    ll[weights <= 0] <- 0
                    if(any(!is.finite(ll))) NaN else sum(weights * ll)
                    ## again: catch extreme cases without warning
                }
            })
        }

        ## gradient (by default) or gradient contributions (sum = FALSE)
        gradfun <- function(par, sum = TRUE, fit = NULL) {
            ## extract fitted means/precisions
            if(is.null(fit)) fit <- fitfun(par, deriv = 1L)
            rval <- with(fit, {
                ## compute gradient contributions
                cbind(
                    phi * (ystar - mustar) * mu.eta(eta) * weights * x,
                    (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
                    phi_mu.eta(phi_eta) * weights * z
                )
            })
            rval[weights <= 0, ] <- 0
            if(sum) colSums(rval) else rval
        }

        ## analytical Hessian (expected information) or covariance matrix (inverse of Hessian)
        hessfun <- function(par, inverse = FALSE, fit = NULL) {
            ## extract fitted means/precisions
            if(is.null(fit)) fit <- fitfun(par, deriv = 2L)
            with(fit, {
                ## auxiliary transformations
                a <- psi1 + psi2
                b <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
                ## compute elements of W
                wbb <- phi^2 * a * mu.eta(eta)^2
                wpp <- b * phi_mu.eta(phi_eta)^2
                wbp <- phi * (mu * a - psi2) * mu.eta(eta) * phi_mu.eta(phi_eta)
                ## compute elements of K
                kbb <- if(k > 0L) crossprod(sqrt(weights) * sqrt(wbb) * x) else crossprod(x)
                kpp <- if(m > 0L) crossprod(sqrt(weights) * sqrt(wpp) * z) else crossprod(z)
                kbp <- if(k > 0L & m > 0L) crossprod(weights * wbp * x, z) else crossprod(x, z)
                ## put together K (= expected information)
                K <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
                if (inverse) chol2inv(chol(K)) else K
                ## previously computed K^(-1) via partitioned matrices, but this appears to be
                ## slower - even for moderately sized problems
                ##   kbb1 <- if(k > 0L) chol2inv(qr.R(qr(sqrt(weights) * sqrt(wbb) * x))) else kbb
                ##   kpp1 <- if(m > 0L) solve(kpp - t(kbp) %*% kbb1 %*% kbp) else kpp
                ##   vcov <- cbind(rbind(kbb1 + kbb1 %*% kbp %*% kpp1 %*% t(kbp) %*% kbb1,
                ##     -kpp1 %*% t(kbp) %*% kbb1), rbind(-kbb1 %*% kbp %*% kpp1, kpp1))
            })
        }

        ## compute biases and adjustment for bias correction/reduction
        biasfun <- function(par, fit = NULL, vcov = NULL) {
            if (is.null(fit)) fit <- fitfun(par, deriv = 2L)
            InfoInv <- if(is.null(vcov)) try(hessfun(par, inverse = TRUE), silent = TRUE) else vcov
            mu <- fit$mu
            with(fit, {
                kappa2 <- psi1 + psi2
                D1 <- mu.eta(eta)
                D2 <- phi_mu.eta(phi_eta)
                D1dash <- d2mu.deta(eta)
                D2dash <- phi_d2mu.deta(phi_eta)
                dpsi1 <-  psigamma(mu * phi, 2)       ## potentially move to fitfun() when we add support for
                dpsi2 <-  psigamma((1 - mu) * phi, 2) ## observed information (as opposed to expected)
                kappa3 <- dpsi1 - dpsi2
                psi3 <- psigamma(phi, 1)
                dpsi3 <- psigamma(phi, 2)
                ## PQsum produces the the adustments to the score functions and is suggested for iteration
                PQsum <- function(t) {
                    if (t <= k)  {
                        Xt <- x[,t]
                        bb <- if (k > 0L)
                                  crossprod(x, weights * phi^2 * D1 * (phi * D1^2 * kappa3 + D1dash * kappa2) * Xt * x)
                              else
                                  crossprod(x)
                        bg <- if ((k > 0L) & (m > 0L))
                                  crossprod(x, weights * phi * D1^2 * D2 * (mu * phi * kappa3 + phi * dpsi2 + kappa2) * Xt * z)
                              else
                                  crossprod(x, z)
                        gg <- if (m > 0L)
                                  crossprod(z, weights * phi * D1 * D2^2 * (mu^2 * kappa3 - dpsi2 + 2 * mu * dpsi2) * Xt * z) +
                                      crossprod(z, weights * phi * D1 * D2dash * (mu * kappa2 - psi2) * Xt * z)
                              else
                                  crossprod(z)
                    } else {
                        Zt <- z[, t - k]
                        bb <- if (k > 0L)
                                  crossprod(x, weights * phi * D2 * (phi * D1^2 * mu * kappa3 + phi * D1^2 * dpsi2 + D1dash * mu * kappa2 - D1dash * psi2) * Zt * x)
                              else
                                  crossprod(x)
                        bg <- if ((k > 0L) & (m > 0L))
                                  crossprod(x, weights * D1 * D2^2 * (phi * mu^2 * kappa3 + phi * (2 * mu - 1) * dpsi2 + mu * kappa2 - psi2) * Zt * z)
                              else
                              crossprod(x, z)
                        gg <- if (m > 0L)
                                  crossprod(z, weights * D2^3 * (mu^3 * kappa3 + (3 * mu^2 - 3 * mu + 1) * dpsi2 - dpsi3) * Zt * z) +
                                      crossprod(z, weights * D2dash * D2 * (mu^2 * kappa2 + (1 - 2 * mu) * psi2 - psi3) * Zt * z)
                              else
                                  crossprod(z)
                    }
                    pq <- rbind(cbind(bb, bg), cbind(t(bg), gg))
                    sum(diag(InfoInv %*% pq))/2
                }
                if (inherits(InfoInv, "try-error")) {
                    bias <- adjustment <- rep.int(NA_real_, k + m)
                }
                else {
                    adjustment <- sapply(1:(k + m), PQsum)
                    bias <- - InfoInv %*% adjustment
                }
                list(bias = bias, adjustment = adjustment)
            })
        }
    }

    ##  if(dist != "beta") {
    else {
        ## distribution d/p functions
        dfun <- if(dist == "xbeta") {
                    function(x, mu, phi, nu, ...) dxbeta(x, mu = mu, phi = phi, nu = nu, ...)
                }
                else {
                    quadrule <- quadtable(nquad = quad)
                    function(x, mu, phi, nu, ...) dxbetax(x, mu = mu, phi = phi, nu = nu, quad = quadrule, ...)
                }

        ## set up (censored) log-likelihood
        loglikfun <- function(par, fit = NULL) {
            ## extract fitted parameters
            if(is.null(fit)) fit <- fitfun(par)
            ## FIXME: Temporary edge case control to go around the
            ## stopifnot's in dxbeta, dxbetax, when simple_formula is
            ## TRUE (link.phi is "identity") and the optimizer tries
            ## negative phi's
            with(fit, {
                if(any(!is.finite(phi)) | any(phi < 0))
                    NaN
                else {
                    rval <- dfun(y, mu = mu, phi = phi, nu = nu, log = TRUE)
                    sum(weights * rval)
                }
            })
        }

        ## FIXME: Numerics are not yet very reliable. Accurate evaluation of h3f2 is still the issue here...
        gradfun_xbeta <- function(par, sum = TRUE, fit = NULL) {
            ## extract fitted means/precisions
            if(is.null(fit)) fit <- fitfun(par, deriv = 3L)
            with(fit, {
                ynu <- (y + nu)/(1 + 2 * nu)
                ystarnu <- qlogis(ynu)
                ## Compute gradient contributions from observations in (0, 1) (beta regression scores)
                grad_l01 <- cbind(
                    phi * (ystarnu - mustar) * mu.eta(eta) * weights * x,
                    (mu * (ystarnu - mustar) + log(1 - ynu) - digamma((1 - mu) * phi) + digamma(phi)) *
                    phi_mu.eta(phi_eta) * weights * z,
                    ((mu * phi - 1)/(y + nu) + ((1 - mu) * phi - 1)/(1 - y + nu) - 2 * (phi - 1)/(1 + 2 * nu)) * weights * nu)
                nu_low <- nu/(1 + 2*nu)
                nu_upp <- (1 + nu)/(1 + 2*nu)
                dlow <- dbeta(nu_low, shape1, shape2)
                plow <- pbeta(nu_low, shape1, shape2)
                dupp <- dbeta(nu_upp, shape1, shape2)
                pupp <- pbeta(nu_upp, shape1, shape2)
                ## Fs1 - Fs4 can explode below and floating point error propagates
                Fs1 <- h3f2(shape1, shape2, nu_low, n, maxiter = 10000, eps = 0)
                Fs2 <- h3f2(shape1, shape2, nu_upp, n, maxiter = 10000, eps = 0)
                Fs3 <- h3f2(shape2, shape1, nu_low, n, maxiter = 10000, eps = 0)
                Fs4 <- h3f2(shape2, shape1, nu_upp, n, maxiter = 10000, eps = 0)
                delta1low <- plow * (d12 - d1 + log(nu_low)) - nu_low^shape1 * Fs1 / (shape1^2 * b12)
                delta2low <- (1 - plow) * (d2 - d12 - log(nu_upp)) + nu_upp^shape2 * Fs4 / (shape2^2 * b12)
                delta1upp <- pupp * (d12 - d1 + log(nu_upp)) - nu_upp^shape1 * Fs2 / (shape1^2 * b12)
                delta2upp <- (1 - pupp) * (d2 - d12 - log(nu_low)) + nu_low^shape2 * Fs3 / (shape2^2 * b12)
                ## Tested the above with numerical derivatives; look ok
                grad_l0 <- cbind(
                    phi * mu.eta(eta) * (delta1low - delta2low) * weights * x / plow,
                    phi_mu.eta(phi_eta) * (delta1low * mu + delta2low * (1 - mu)) * weights * z / plow,
                    dlow/(plow * (1 + 2 * nu)^2) * weights * nu ## case weights here
                )
                grad_l1 <- cbind(
                    phi * mu.eta(eta) * (delta2upp - delta1upp) * weights * x / (1 - pupp),
                    - phi_mu.eta(phi_eta) * (delta1upp * mu + delta2upp * (1 - mu)) * weights * z / (1 - pupp),
                    dupp/((1 - pupp) * (1 + 2 * nu)^2) * weights * nu ## case weights here
                )
                grad_l0[!indices0, ] <- 0
                grad_l1[!indices1, ] <- 0
                grad_l01[!indices01, ] <- 0
                out <- grad_l0 + grad_l1 + grad_l01
                out <- if (estnu) out else out[, 1:(k + m)]
                if (sum) {
                    colSums(out)
                }
                else {
                    out
                }
            })
        }

        if (dist == "xbeta") {
            gradfun <- gradfun_xbeta
        }
        if (dist == "xbetax") {
            gradfun <- function(par, sum = TRUE, fit = NULL) {
                if (is.null(fit)) {
                    fit <- fitfun(par, deriv = 3L)
                }
                with(fit, {
                    dens <- apply(quadrule, 1, function(rule) {
                        e <- rule[1] * nu
                        rule[2] * dxbeta(y, mu, phi, nu = e, log = FALSE)
                    })
                    tdens <- rowSums(dens)
                    obsders <- lapply(seq.int(quad), function(inds) {
                        current_nu <- fit$nu <- quadrule[inds, 1]*nu
                        par[k + m + 1] <- log(current_nu)
                        out <- gradfun_xbeta(par, fit = fit, sum = FALSE)
                        dens[, inds]*out/tdens
                    })
                    out <- Reduce("+", obsders)
                    if (sum) {
                        colSums(out)
                    }
                    else {
                        out
                    }
                })
            }
        }
    }

    ## optimize likelihood
    if (method == "nlminb") {
        stopifnot(requireNamespace("numDeriv"))
        if("maxit" %in% control) {
          if(is.null(control$iter.max)) control$iter.max <- control$maxit
          control$maxit <- NULL
        }
        if("reltol" %in% control) {
          if(is.null(control$rel.tol)) control$rel.tol <- control$reltol
          control$reltol <- NULL
        }
        opt <- nlminb(start = start, objective = function(par, ...) -loglikfun(par, ...),
                      gradient = if (gradient) function(par, ...) -gradfun(par, ...) else NULL,
                      control = control)
        opt$hessian <- numDeriv::hessian(loglikfun, opt$par)
    }
    else {
        opt <- optim(par = start, fn = loglikfun, gr = if (gradient) gradfun else NULL,
                     method = method, hessian = hessian, control = control)
    }
    par <- opt$par

    ## conduct further (quasi) Fisher scoring to move ML derivatives
    ## even further to zero or conduct bias reduction
    ## (suppressed if fsmaxit = 0 or if only numerical optim result desired)
    if(type == "BR" & fsmaxit <= 0) warning("BR cannot be performed with fsmaxit <= 0")
    step <- .Machine$integer.max
    iter <- 0
    if(fsmaxit > 0 & !(hessian & type == "ML"))
    {
        for (iter in 1:fsmaxit) {
            stepPrev <- step
            stepFactor <- 0
            testhalf <- TRUE
            while (testhalf & stepFactor < 11) {
                fit <- fitfun(par, deriv = 2L)
                scores <- gradfun(par, fit = fit)
                InfoInv <- try(hessfun(par, fit = fit, inverse = TRUE))
                if(failedInv <- inherits(InfoInv, "try-error")) {
                    warning("failed to invert the information matrix: iteration stopped prematurely")
                    break
                }
                bias <- if(type == "BR") biasfun(par, fit = fit, vcov = InfoInv)$bias else 0
                par <- par + 2^(-stepFactor) * (step <- InfoInv %*% scores - bias)
                stepFactor <- stepFactor + 1
                testhalf <- drop(crossprod(stepPrev) < crossprod(step))
            }
            if (failedInv | (all(abs(step) < fstol))) {
                break
            }
        }
    }

    ## check whether both optim() and manual iteration converged IK:
    ## modified the condition a bit... optim might fail to converge but
    ## if additional iteration are requested Fisher scoring might get
    ## there
    if((fsmaxit == 0 & opt$convergence > 0) | (iter >= fsmaxit & fsmaxit > 0)) {
        converged <- FALSE
        warning("optimization failed to converge")
    } else {
        converged <- TRUE
    }

    ## conduct single bias correction (if BC selected) else do not
    ## estimate the first order biases
    if(type == "BC") {
        bias <- as.vector(biasfun(par)$bias)
        par <- par - bias
    }
    else {
        bias <- rep.int(NA_real_, k + m + estnu)
    }

    ## extract fitted values/parameters
    fit <- fitfun(par, deriv = 3L)
    beta <- fit$beta
    gamma <- fit$gamma
    eta <- fit$eta
    mu <- fit$mu
    phi <- fit$phi
    nu <- if(!estnu) nu else as.vector(exp(par[k + m + 1]))

    ## log-likelihood/gradients/covariance matrix at optimized parameters
    ll <- loglikfun(par, fit = fit)


    ## No need to evaluate ef below.
    if (gradient) {
        ef <- gradfun(par, fit = fit, sum = FALSE)
    }
    else {
        stopifnot(requireNamespace("numDeriv"))
        ef <- numDeriv::grad(loglikfun, par)
    }

    vcov <- if (hessian && (type == "ML")) solve(-as.matrix(opt$hessian)) else hessfun(fit = fit, inverse = TRUE)

    ## R-squared
    wcor <- function(x, y, weights = NULL) {
      if(is.null(weights) || identical(rep.int(1, length(x)), weights)) return(cor(x, y))
      x <- x[weights > 0]
      y <- y[weights > 0]
      w <- weights[weights > 0]/sum(weights)
      x <- x - sum(x * w)
      x <- x/sqrt(sum(w * x^2))
      y <- y - sum(y * w)
      y <- y/sqrt(sum(w * y^2))
      sum(w * x * y)
    }
    pseudor2 <- if(dist != "beta" || var(eta[weights > 0]) * var(ystar[weights > 0]) <= 0) NA else wcor(eta, linkfun(y), weights)^2

    ## names
    names(beta) <- colnames(x)
    names(gamma) <- if(phi_const & phi_linkstr == "identity") "(phi)" else colnames(z)
    rownames(vcov) <- colnames(vcov) <- names(bias) <- c(colnames(x),
                                                         if(phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"),
                                                         if(estnu) "Log(nu)" else NULL)

    marg_e <- switch(dist,
                     "beta" = mu,
                     "xbeta" = vapply(seq.int(n), function(i) mean_xbeta(mu[i], phi[i], nu), 0.0),
                     "xbetax" = vapply(seq.int(n), function(i) mean_xbetax(mu[i], phi[i], nu, quad), 0.0))

    ## set up return value
    rval <- list(
        coefficients = list(mu = beta, phi = gamma),
        residuals = y - marg_e,
        fitted.values = structure(marg_e, .Names = names(y)),
        type = type,
        dist = dist,
        optim = opt,
        method = method,
        control = ocontrol,
        scoring = iter,
        start = start,
        weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
        offset = list(mu = if(identical(offset[[1L]], rep.int(0, n))) NULL else offset[[1L]],
                      phi = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
        n = n,
        nobs = nobs,
        df.null = nobs - 2 - estnu,
        df.residual = nobs - k - m - estnu,
        phi = phi_full,
        nu = nu,
        loglik = ll,
        vcov = vcov,
        bias = bias,
        pseudo.r.squared = pseudor2,
        link = list(mu = linkobj, phi = phi_linkobj),
        converged = converged,
        grad = ef
    )
    if(estnu) rval$coefficients$nu <- c("Log(nu)" = log(nu))
    if(dist == "beta") {
        for(n in intersect(names(rval), fix_names_mu_phi)) names(rval[[n]])[1L:2L] <- c("mean", "precision")
    }
    return(rval)
}

print.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## unify list component names
    if(is.null(x$dist) || (x$dist == "beta")) {
        for(n in intersect(names(x), fix_names_mu_phi)) names(x[[n]])[1L:2L] <- c("mu", "phi")
        mp <- c("mean", "precision")
    } else {
        mp <- c("mu", "phi")
    }

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

    if(!x$converged) {
        cat("model did not converge\n")
    } else {
        if(length(x$coefficients$mu)) {
            cat(sprintf("Coefficients (%s model with %s link):\n", mp[1L], x$link$mu$name))
            print.default(format(x$coefficients$mu, digits = digits), print.gap = 2, quote = FALSE)
            cat("\n")
        } else cat(sprintf("No coefficients (in %s model)\n\n", mp[1L]))
        if(x$phi) {
            if(length(x$coefficients$phi)) {
                cat(sprintf("Phi coefficients (%s model with %s link):\n", mp[2L], x$link$phi$name))
                print.default(format(x$coefficients$phi, digits = digits), print.gap = 2, quote = FALSE)
                cat("\n")
            } else cat(sprintf("No coefficients (in %s model)\n\n", mp[2L]))
        }
    }
    if(!is.null(x$dist) && (x$dist != "beta")) {
        cat(sprintf("Exceedence parameter (extended-support %s model)\nnu: %s\n\n",
                    x$dist,
                    round(x$nu, digits = digits)))
    }

    invisible(x)
}

summary.betareg <- function(object, phi = NULL, type = "quantile", ...)
{
    ## unify list component names
    if(is.null(object$dist) || (object$dist == "beta")) {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mu", "phi")
    }

    ## treat phi as full model parameter?
    if(!is.null(phi)) object$phi <- phi

    ## residuals
    type <- match.arg(type, c("quantile", "pearson", "deviance", "response", "weighted", "sweighted", "sweighted2"))
    object$residuals <- residuals(object, type = type)
    object$residuals.type <- type

    ## extend coefficient table
    k <- length(object$coefficients$mu)
    m <- length(object$coefficients$phi)
    cf <- as.vector(do.call("c", object$coefficients))
    se <- sqrt(diag(object$vcov))
    cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
    colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    nu <- if("nu" %in% names(object$coefficients)) cf[m + k + 1, , drop = FALSE] else NULL
    cf <- list(
        mu = cf[seq.int(length.out = k), , drop = FALSE],
        phi = cf[seq.int(length.out = m) + k, , drop = FALSE]
    )
    rownames(cf$mu) <- names(object$coefficients$mu)
    rownames(cf$phi) <- names(object$coefficients$phi)
    if(!is.null(nu)) {
        cf$nu <- nu
        rownames(cf$nu) <- names(object$coefficients$nu)
    }
    object$coefficients <- cf

    ## number of iterations
    mytail <- function(x) x[length(x)]
    if (object$method == "nlminb") {
        object$iterations <- c("nlminb" = as.vector(object$optim$iterations), "scoring" = as.vector(object$scoring))
    } else {
        object$iterations <- c("optim" = as.vector(mytail(na.omit(object$optim$counts))), "scoring" = as.vector(object$scoring))
    }


    ## delete some slots
    object$fitted.values <- object$terms <- object$model <- object$y <-
        object$x <- object$levels <- object$contrasts <- object$start <- NULL

    ## restore old list component names for backward compatibility
    if(is.null(object$dist) || (object$dist == "beta")) {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mean", "precision")
    }

    ## return
    class(object) <- "summary.betareg"
    object
}

print.summary.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## unify list component names
    if(is.null(x$dist) || (x$dist == "beta")) {
        for(n in intersect(names(x), fix_names_mu_phi)) names(x[[n]])[1L:2L] <- c("mu", "phi")
        mp <- c("mean", "precision")
    } else {
        mp <- c("mu", "phi")
    }

    cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

    if(!x$converged) {
        cat("model did not converge\n")
    } else {
        types <- c("quantile", "pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
        Types <- c("Quantile residuals", "Pearson residuals", "Deviance residuals", "Raw response residuals",
                   "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")
        if(!is.null(x$dist) && (x$dist != "beta")) Types[1L] <- "Randomized quantile residuals"
        cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
        print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
                        .Names = c("Min", "1Q", "Median", "3Q", "Max")))

        if(NROW(x$coefficients$mu)) {
            cat(sprintf("\nCoefficients (%s model with %s link):\n", mp[1L], x$link$mu$name))
            printCoefmat(x$coefficients$mu, digits = digits, signif.legend = FALSE)
        } else cat("\nNo coefficients (in mean model)\n")

        if(x$phi) {
            if(NROW(x$coefficients$phi)) {
                cat(sprintf("\nPhi coefficients (%s model with %s link):\n", mp[2L], x$link$phi$name))
                printCoefmat(x$coefficients$phi, digits = digits, signif.legend = FALSE)
            } else cat("\nNo coefficients (in precision model)\n")
        }

        if(!is.null(x$coefficients$nu)) {
            cat("\nExceedence parameter (extended-support xbetax model):\n")
            printCoefmat(x$coefficients$nu, digits = digits, signif.legend = FALSE)
        }

        if(getOption("show.signif.stars") && any(do.call("rbind", x$coefficients)[, 4L] < 0.1, na.rm = TRUE))
            cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")


        if(!is.null(x$dist) && (x$dist != "beta")) {
            cat(sprintf("\nExceedence parameter nu: %s", round(x$nu, digits = digits)))
        }
        cat("\nType of estimator:", x$type, switch(x$type,
                                                   "ML" = "(maximum likelihood)",
                                                   "BC" = "(bias-corrected)",
                                                   "BR" = "(bias-reduced)"))
        cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
            "on", sum(sapply(x$coefficients, NROW)), "Df")
        if(!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
        if(x$iterations[2L] > 0) {
            scoring_type <- switch(x$type,
                                   "ML" = "(Fisher scoring)",
                                   "BR" = "(quasi Fisher scoring)",
                                   "BC" = "(Fisher scoring)")
            cat(paste("\nNumber of iterations:", x$iterations[1L],
                      sprintf("(%s) +", x$method), x$iterations[2L], paste(scoring_type, "\n")))
        } else {
            cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations[1L], "\n"))
        }
    }

    invisible(x)
}

predict.betareg <- function(object, newdata = NULL,
                            type = c("response", "link", "precision", "variance", "parameters", "distribution", "density", "probability", "quantile"),
                            na.action = na.pass, at = 0.5, elementwise = NULL, ...)
{
    ## types of predictions
    type <- match.arg(type[1L], c(
        "response", "mean",
        "link",
        "precision",
        "variance",
        "density", "pdf",
        "probability", "cdf",
        "quantile",
        "distribution",
        "parameters"))
    if(type == "mean") type <- "response"
    if(type == "cdf") type <- "probability"
    if(type == "pdf") type <- "density"

    ## for distribution require distributions3 infrastructure
    if(type == "distribution") stopifnot(requireNamespace("distributions3"))

    ## unify list component names
    if(is.null(object$dist)) object$dist <- "beta"
    if(object$dist == "beta") {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mu", "phi")
    }
    if(is.null(object$nu)) object$nu <- NA_real_

    ## set up function that computes prediction from model parameters
    linkfun <- object$link$mu$linkfun
    fun <- if(object$dist == "beta") {
               switch(type,
                      "response" = function(pars) pars$mu,
                      "link" = function(pars) linkfun(pars$mu),
                      "precision" = function(pars) pars$phi,
                      "variance" = function(pars) pars$mu * (1 - pars$mu)/(1 + pars$phi),
                      "parameters" = function(pars) {pars$nu <- NULL; pars},
                      "distribution" = function(pars) BetaR(mu = pars$mu, phi = pars$phi),
                      "density" = function(x, pars, ...) dbetar(x, mu = pars$mu, phi = pars$phi, ...),
                      "probability" = function(q, pars, ...) pbetar(q, mu = pars$mu, phi = pars$phi, ...),
                      "quantile" = function(p, pars, ...) qbetar(p, mu = pars$mu, phi = pars$phi, ...)
                      )
           } else if(object$dist == "xbeta") {
               switch(type,
                      "response" = function(pars) apply(pars, 1, function(x) mean_xbeta(x["mu"], x["phi"], x["nu"])),
                      "link" = function(pars) linkfun(pars$mu),
                      "precision" = function(pars) pars$phi,
                      "variance" = function(pars) apply(pars, 1, function(x) var_xbeta(x["mu"], x["phi"], x["nu"])),
                      "parameters" = function(pars) pars,
                      "distribution" = function(pars) XBeta(mu = pars$mu, phi = pars$phi, nu = pars$nu),
                      "density" = function(x, pars, ...) dxbeta(x, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...),
                      "probability" = function(q, pars, ...) pxbeta(q, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...),
                      "quantile" = function(p, pars, ...) qxbeta(p, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...)
                      )
           } else if(object$dist == "xbetax") {
               switch(type,
                      "response" = function(pars) apply(pars, 1, function(x) mean_xbetax(x["mu"], x["phi"], x["nu"], object$control$quad)),
                      "link" = function(pars) linkfun(pars$mu),
                      "precision" = function(pars) pars$phi,
                      "variance" = function(pars) apply(pars, 1, function(x) var_xbetax(x["mu"], x["phi"], x["nu"], object$control$quad)),
                      "parameters" = function(pars) pars,
                      "distribution" = function(pars) XBetaX(mu = pars$mu, phi = pars$phi, nu = pars$nu),
                      "density" = function(x, pars, ...) dxbetax(x, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...),
                      "probability" = function(q, pars, ...) pxbetax(q, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...),
                      "quantile" = function(p, pars, ...) qxbetax(p, mu = pars$mu, phi = pars$phi, nu = pars$nu, ...)
                      )
           }

    if(missing(newdata) || is.null(newdata)) {
        pars <- data.frame(mu = object$fitted.values, phi = NA_real_, nu = object$nu)
        if(!(type %in% c("link"))) {
            if (object$dist != "beta") {
                ## Use the correct mu values
                beta <- object$coefficients$mu
                x <- if(is.null(object$x)) model.matrix(object, model = "mu") else object$x$mu
                offset <- if(is.null(object$offset$mu)) rep.int(0, NROW(x)) else object$offset$mu
                pars$mu <- object$link$mu$linkinv(drop(x %*% beta + offset))
            }
            gamma <- object$coefficients$phi
            z <- if(is.null(object$x)) model.matrix(object, model = "phi") else object$x$phi
            offset <- if(is.null(object$offset$phi)) rep.int(0, NROW(z)) else object$offset$phi
            pars$phi <- object$link$phi$linkinv(drop(z %*% gamma + offset))
        }
    } else {
        tnam <- switch(type,
                       "response" = if(object$dist == "beta") "mu" else "full",
                       "link" = "mu",
                       "precision" = "phi",
                       "full")

        mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
        newdata <- newdata[rownames(mf), , drop = FALSE]
        offset <- list(mu = rep.int(0, nrow(mf)), phi = rep.int(0, nrow(mf)))

        pars <- data.frame(mu = rep.int(NA_real_, nrow(mf)), phi = NA_real_, nu = object$nu)
        rownames(pars) <- rownames(mf)

        if(type != "precision") {
            X <- model.matrix(delete.response(object$terms$mu), mf, contrasts.arg = object$contrasts$mu)
            if(!is.null(object$call$offset)) offset[[1L]] <- offset[[1L]] + eval(object$call$offset, newdata)
            if(!is.null(off.num <- attr(object$terms$mu, "offset"))) {
                for(j in off.num) offset[[1L]] <- offset[[1L]] + eval(attr(object$terms$mu, "variables")[[j + 1L]], newdata)
            }
            pars$mu <- object$link$mu$linkinv(drop(X %*% object$coefficients$mu + offset[[1L]]))
        }
        if(!((object$dist == "beta" && type == "response") || (type == "link"))) {
            Z <- model.matrix(object$terms$phi, mf, contrasts.arg = object$contrasts$phi)
            if(!is.null(off.num <- attr(object$terms$phi, "offset"))) {
                for(j in off.num) offset[[2L]] <- offset[[2L]] + eval(attr(object$terms$phi, "variables")[[j + 1L]], newdata)
            }
            pars$phi <- object$link$phi$linkinv(drop(Z %*% object$coefficients$phi + offset[[2L]]))
        }
    }

    if(type %in% c("response", "link", "precision", "variance", "parameters", "distribution")) {
        ## prediction is just a transformation of the parameters
        rval <- fun(pars, ...)
        if(is.null(dim(rval))) {
          names(rval) <- rownames(pars)
        } else {
          rownames(rval) <- rownames(pars)
        }
    } else {
        ## prediction requires a function that suitably expands 'at'
        ## and then evaluates fun() with predicted parameters as default
        FUN <- function(at, mu = pars$mu, phi = pars$phi, nu = pars$nu, elementwise = NULL, ...) {
            n <- length(mu)
            if(is.null(elementwise)) {
                elementwise <- !( all(length(at) != c(1L, n)) || (is.matrix(at) && NROW(at) == 1L) )
            }
            if(elementwise) {
                if(length(at) == 1L) at <- rep.int(as.vector(at), n)
                if(length(at) != n) stop("for 'elementwise = TRUE' the argument 'at' must either have length 1 or the same as the number of observations in 'newdata'")
                rv <- fun(at, pars = data.frame(mu = mu, phi = phi, nu = nu))
                names(rv) <- rownames(pars)
            } else {
                at <- matrix(rep(at, each = n), nrow = n)
                rv <- fun(as.vector(at), pars = data.frame(mu = rep.int(mu, ncol(at)),
                                                           phi = rep.int(phi, ncol(at)), nu = rep.int(nu, ncol(at))), ...)
                rv <- matrix(rv, nrow = n)
                rownames(rv) <- rownames(pars)
                colnames(rv) <- paste(substr(type, 1L, 1L),
                                      round(at[1L, ], digits = pmax(3L, getOption("digits") - 3L)), sep = "_")
            }
            return(rv)
        }
        rval <- FUN(at, elementwise = elementwise, ...)
    }

    return(rval)
}

coef.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...)
{
    ## unify list component names
    if(is.null(object$dist) || (object$dist == "beta")) {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mu", "phi")
    }

    cf <- object$coefficients

    model <- if(is.null(phi)) {
                 if(missing(model)) ifelse(object$phi, "full", "mu") else fix_model_mu_phi(model)[1L]
             } else {
                 if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
                 ifelse(phi, "full", "mu")
             }

    switch(model,
           "mu" = {
               cf$mu
           },
           "phi" = {
               cf$phi
           },
           "nu" = {
               cf$nu
           },
           "full" = {
               nam1 <- names(cf$mu)
               nam2 <- names(cf$phi)
               nam3 <- names(cf$nu)
               cf <- c(cf$mu, cf$phi, cf$nu)
               names(cf) <- c(nam1, if(identical(nam2, "(phi)")) "(phi)" else paste("(phi)", nam2, sep = "_"), nam3)
               cf
           }
           )
}

vcov.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...)
{
    ## unify list component names
    if(is.null(object$dist) || (object$dist == "beta")) {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mu", "phi")
    }

    vc <- object$vcov
    k <- length(object$coefficients$mu)
    m <- length(object$coefficients$phi)

    model <- if(is.null(phi)) {
                 if(missing(model)) ifelse(object$phi, "full", "mu") else fix_model_mu_phi(model)[1L]
             } else {
                 if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
                 ifelse(phi, "full", "mu")
             }

    switch(model,
           "mu" = {
               vc[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
           },
           "phi" = {
               vc <- vc[seq.int(length.out = m) + k, seq.int(length.out = m) + k, drop = FALSE]
               colnames(vc) <- rownames(vc) <- names(object$coefficients$phi)
               vc
           },
           "nu" = {
               vc[m + k + 1, m + k + 1, drop = FALSE]
           },
           "full" = {
               vc
           }
           )
}

bread.betareg <- function(x, phi = NULL, ...) {
    vcov(x, phi = phi) * x$nobs
}

estfun.betareg <- function(x, phi = NULL, ...)
{
    ## unify list component names
    if(is.null(x$dist) || (x$dist == "beta")) {
        for(n in intersect(names(x), fix_names_mu_phi)) names(x[[n]])[1L:2L] <- c("mu", "phi")
    } else {
        stop("not yet implemented")
    }

    ## extract response y and regressors X and Z
    y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
    xmat <- if(is.null(x$x)) model.matrix(x, model = "mu") else x$x$mu
    zmat <- if(is.null(x$x)) model.matrix(x, model = "phi") else x$x$phi
    offset <- x$offset
    if(is.null(offset[[1L]])) offset[[1L]] <- rep.int(0, NROW(xmat))
    if(is.null(offset[[2L]])) offset[[2L]] <- rep.int(0, NROW(zmat))
    wts <- weights(x)
    if(is.null(wts)) wts <- 1
    phi_full <- if(is.null(phi)) x$phi else phi

    ## extract coefficients
    beta <- x$coefficients$mu
    gamma <- x$coefficients$phi

    ## compute y*
    ystar <- qlogis(y)

    ## compute mu*
    eta <- as.vector(xmat %*% beta + offset[[1L]])
    phi_eta <- as.vector(zmat %*% gamma + offset[[2L]])
    mu <- x$link$mu$linkinv(eta)
    phi <- x$link$phi$linkinv(phi_eta)
    mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

    ## compute scores of beta
    rval <- phi * (ystar - mustar) * as.vector(x$link$mu$mu.eta(eta)) * wts * xmat

    ## combine with scores of phi
    if(phi_full) {
        rval <- cbind(rval,
        (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
        as.vector(x$link$phi$mu.eta(phi_eta)) * wts * zmat)
        colnames(rval) <- names(coef(x, phi = phi_full))
    }
    attr(rval, "assign") <- NULL

    ##
    if(x$type == "BR") {
        nobs <- nrow(xmat)
        k <- ncol(xmat)
        m <- ncol(zmat)
        InfoInv <- x$vcov
        D1 <- x$link$mu$mu.eta(eta)
        D2 <- x$link$phi$mu.eta(phi_eta)
        D1dash <- x$link$mu$d2mu.deta(eta)
        D2dash <- x$link$phi$d2mu.deta(phi_eta)
        psi2 <- psigamma((1 - mu) * phi, 1)
        dpsi1 <-  psigamma(mu * phi, 2)
        dpsi2 <-  psigamma((1 - mu) * phi, 2)
        kappa2 <- psigamma(mu * phi, 1) + psi2
        kappa3 <- dpsi1 - dpsi2
        psi3 <- psigamma(phi, 1)
        dpsi3 <- psigamma(phi, 2)
        PQ <- function(t) {
            prodfun <- function(mat1, mat2) {
                sapply(seq_len(nobs), function(i) tcrossprod(mat1[i,], mat2[i,]), simplify = FALSE)
            }
            if (t <= k)  {
                Xt <- xmat[,t]
                bb <- if (k > 0L) {
                          bbComp <- wts * phi^2 * D1 * (phi * D1^2 * kappa3 + D1dash * kappa2) * Xt * xmat
                          prodfun(xmat, bbComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, k, k))
                bg <- if ((k > 0L) & (m > 0L)) {
                          bgComp <- wts * phi * D1^2 * D2 * (mu * phi * kappa3 + phi * dpsi2 + kappa2) * Xt * zmat
                          prodfun(xmat, bgComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, k, m))
                gg <- if (m > 0L) {
                          ggComp <- wts * phi * D1 * D2^2 * (mu^2 * kappa3 - dpsi2 + 2 * mu * dpsi2) * Xt * zmat +
                              wts * phi * D1 * D2dash * (mu * kappa2 - psi2) * Xt * zmat
                          prodfun(zmat, ggComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, m, m))
            } else {
                Zt <- zmat[, t - k]
                bb <- if (k > 0L) {
                          bbComp <- wts * phi * D2 * (phi * D1^2 * mu * kappa3 + phi * D1^2 * dpsi2 + D1dash * mu * kappa2 - D1dash * psi2) * Zt * xmat
                          prodfun(xmat, bbComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, k, k))
                bg <- if ((k > 0L) & (m > 0L)) {
                          bgComp <- wts * D1 * D2^2 * (phi * mu^2 * kappa3 + phi * (2 * mu - 1) * dpsi2 + mu * kappa2 - psi2) * Zt * zmat
                          prodfun(xmat, bgComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, k, m))
                gg <- if (m > 0L) {
                          ggComp <- wts * D2^3 * (mu^3 * kappa3 + (3 * mu^2 - 3 * mu + 1) * dpsi2 - dpsi3) * Zt * zmat +
                              wts * D2dash * D2 * (mu^2 * kappa2 + (1 - 2 * mu) * psi2 - psi3) * Zt * zmat
                          prodfun(zmat, ggComp)
                      }
                      else
                          sapply(1:nobs, function(x) matrix(0, m, m))
            }
            sapply(seq_len(nobs), function(i)
                sum(diag(InfoInv %*% rbind(cbind(bb[[i]], bg[[i]]), cbind(t(bg[[i]]), gg[[i]]))))/2,
                simplify = TRUE)
        }
        if (inherits(InfoInv, "try-error")) {
            adjustment <- rep.int(NA_real_, k + m)
        }
        else
            adjustment <- sapply(1:(k + m), PQ)
        rval <- rval + adjustment
    }
    return(rval)
}

coeftest.betareg <- function(x, vcov. = NULL, df = Inf, ...)
    coeftest.default(x, vcov. = vcov., df = df, ...)

logLik.betareg <- function(object, ...) {
    structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}

terms.betareg <- function(x, model = c("mean", "precision"), ...)
{
    names(x$terms)[1L:2L] <- c("mu", "phi")
    x$terms[[fix_model_mu_phi(model)[1L]]]
}

model.frame.betareg <- function(formula, ...) {
    if(!is.null(formula$model)) return(formula$model)
    formula$terms <- formula$terms$full
    formula$call$formula <- formula$formula <- formula(formula$terms)
    NextMethod()
}

model.matrix.betareg <- function(object, model = c("mean", "precision"), ...) {
    model <- fix_model_mu_phi(model)[1L]
    for(n in names(object)[names(object) %in% c("x", "terms", "contrasts")]) names(object[[n]])[1L:2L] <- c("mu", "phi")
    rval <- if(!is.null(object$x[[model]])) object$x[[model]]
            else model.matrix(object$terms[[model]], model.frame(object), contrasts.arg = object$contrasts[[model]])
    return(rval)
}

residuals.betareg <- function(object,
                              type = c("quantile", "deviance", "pearson", "response", "weighted", "sweighted", "sweighted2"), ...)
{
    ## unify list component names
    type <- match.arg(type)
    if(is.null(object$dist)) object$dist <- "beta"
    if(object$dist == "beta") {
        for(n in intersect(names(object), fix_names_mu_phi)) names(object[[n]])[1L:2L] <- c("mu", "phi")
    } else {
        ## FIXME: Support more types of residuals, at least "deviance"?
        if (!(type %in% c("quantile", "response", "pearson")))
            stop(sprintf("only types 'quantile', 'response', and 'pearson' are implemented for '%s'", object$dist))
    }

    ## raw response residuals and desired type
    res <- object$residuals
    if(type == "response") return(res)

    ## extract fitted information
    y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
    wts <- weights(object)
    if(is.null(wts)) wts <- 1
    pars <- predict(object, type = "parameters")
    mu <- pars$mu
    phi <- pars$phi
    nu <- pars$nu

    res <- switch(type,

                  "pearson" = {
                      margvar <- switch(object$dist,
                        "beta"   = mu * (1 - mu)/(1 + phi),
                        "xbeta"  = vapply(seq_along(mu), function(i) var_xbeta(mu[i], phi[i], nu[i]), 0.0),
                        "xbetax" = vapply(seq_along(mu), function(i) var_xbetax(mu[i], phi[i], nu[i], object$control$quad), 0.0)
                      )
                      sqrt(wts) * res / sqrt(margvar)
                  },
                  
                  "quantile" = {
                    ## probability integral transform
                    pit <- switch(object$dist,
                      "beta"   =  pbetar(y, mu = mu, phi = phi),
                      "xbeta"  =  pxbeta(y, mu = mu, phi = phi, nu = nu),
                      "xbetax" = pxbetax(y, mu = mu, phi = phi, nu = nu, quad = object$control$quad)
                    )
                    ## boundary observations?
                    if(any(i0 <- y <= 0)) pit[i0] <- pit[i0] * runif(sum(i0))
                    if(any(i1 <- y >= 1)) {
                      pit1 <- switch(object$dist,
                        "beta"   =  pbetar(y, mu = mu, phi = phi, lower.tail = FALSE),
                        "xbeta"  =  pxbeta(y, mu = mu, phi = phi, nu = nu, lower.tail = FALSE),
                        "xbetax" = pxbetax(y, mu = mu, phi = phi, nu = nu, lower.tail = FALSE, quad = object$control$quad)
                      )
                      pit[i1] <- 1 - (pit1[i1] * runif(sum(i1)))
                    }
                    sqrt(wts) * qnorm(pit)
                  },

                  "deviance" = {
                      ll <- function(mu, phi)
                      (lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
                       (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y))
                      sqrt(wts) * sign(res) * sqrt(2 * abs(ll(y, phi) - ll(mu, phi)))
                  },

                  "weighted" = {
                      ystar <- qlogis(y)
                      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
                      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
                      sqrt(wts) * (ystar - mustar) / sqrt(phi * v)
                  },

                  "sweighted" = {
                      ystar <- qlogis(y)
                      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
                      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
                      sqrt(wts) * (ystar - mustar) / sqrt(v)
                  },

                  "sweighted2" = {
                      ystar <- qlogis(y)
                      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
                      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
                      sqrt(wts) * (ystar - mustar) / sqrt(v * (1 - hatvalues(object)))
                  })

    return(res)
}

cooks.distance.betareg <- function(model, ...)
{
    ## unify list component names
    if(is.null(model$dist) || (model$dist == "beta")) {
        for(n in intersect(names(model), fix_names_mu_phi)) names(model[[n]])[1L:2L] <- c("mu", "phi")
    } else {
        stop("not yet implemented for extended-support beta regression")
    }

    h <- hatvalues(model)
    k <- length(model$coefficients$mu)
    res <- residuals(model, type = "pearson")
    h * (res^2)/(k * (1 - h)^2)
}

update.betareg <- function (object, formula., ..., evaluate = TRUE)
{
    call <- object$call
    if(is.null(call)) stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
    if(length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if(any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if(evaluate) eval(call, parent.frame())
    else call
}

simulate.betareg <- function(object, nsim = 1, seed = NULL, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    p <- predict(object, type = "parameter")
    n <- nrow(p)
    nm <- rownames(p)
    if(is.null(object$dist)) object$dist <- "beta"
    s <- switch(object$dist,    
      "beta"   = replicate(nsim, rbetar (n, mu = p$mu, phi = p$phi)),
      "xbeta"  = replicate(nsim, rxbeta (n, mu = p$mu, phi = p$phi, nu = p$nu)),
      "xbetax" = replicate(nsim, rxbetax(n, mu = p$mu, phi = p$phi, nu = p$nu)))
    s <- as.data.frame(s)
    names(s) <- paste("sim", seq_len(nsim), sep = "_")
    if (!is.null(nm)) row.names(s) <- nm
    attr(s, "seed") <- RNGstate
    return(s)
}

prodist.betareg <- function(object, newdata = NULL, na.action = na.pass, ...) {
  if(is.null(object$dist)) object$dist <- "beta"
  predict(object, newdata = newdata, na.action = na.action, type = "distribution", ...)
}
