betareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"),
                    link.phi = NULL,
 		    control = betareg.control(...),
		    model = TRUE, y = TRUE, x = FALSE, ...)
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
  if(length(formula)[2] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  mtZ <- delete.response(terms(formula, data = data, rhs = 2))
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(!(min(Y) > 0 & max(Y) < 1)) stop("invalid dependent variable, all observations must be in (0, 1)")

  ## convenience variables
  n <- length(Y)

  ## links
  if(is.character(link)) link <- match.arg(link)
  if(is.null(link.phi)) link.phi <- if(simple_formula) "identity" else "log"
  if(is.character(link.phi)) link.phi <- match.arg(link.phi, c("identity", "log", "sqrt"))

  ## weights
  weights <- model.weights(mf)
  if(is.null(weights)) weights <- 1
  if(length(weights) == 1) weights <- rep(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  
  ## offset
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if(length(offset) == 1) offset <- rep(offset, n)
  offset <- as.vector(offset)

  ## call the actual workhorse: betareg.fit()
  rval <- betareg.fit(X, Y, Z, weights, offset, link, link.phi, control)

  ## further model information
  rval$call <- cl
  rval$formula <- oformula
  rval$terms <- list(mean = mtX, precision = mtZ, full = mt)
  rval$levels <- list(mean = .getXlevels(mtX, mf), precision = .getXlevels(mtZ, mf), full = .getXlevels(mt, mf))
  rval$contrasts <- list(mean = attr(X, "contrasts"), precision = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(mean = X, precision = Z)

  class(rval) <- "betareg"  
  return(rval)
}

betareg.control <- function(phi = TRUE,
  method = "BFGS", maxit = 5000, hessian = FALSE, trace = FALSE, start = NULL, ...)
{
  rval <- list(phi = phi, method = method, maxit = maxit, hessian = hessian, trace = trace, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.2)
  rval
}

betareg.fit <- function(x, y, z = NULL, weights = NULL, offset = NULL,
  link = "logit", link.phi = "log", control = betareg.control())
{
  ## response and regressor matrix
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep(1, n)
  if(is.null(offset)) offset <- rep(0, n)
  if(is.null(z)) {
    m <- 1
    z <- matrix(1, ncol = m, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
    phi_const <- TRUE
  } else {
    m <- NCOL(z)
    phi_const <- (m == 1L) && isTRUE(all.equal(z[,1], rep(1, n)))
  }

  ## link processing
  if(is.character(link)) {
    linkstr <- link
    linkobj <- if(linkstr != "loglog") make.link(linkstr) else {
      structure(list(
        linkfun = function(mu) -log(-log(mu)),
        linkinv = function(eta) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
        mu.eta = function(eta) {
          eta <- pmin(eta, 700)
          pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
        },
        valideta = function(eta) TRUE,
        name = "loglog"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link
    linkstr <- link$name
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta
  phi_linkstr <- link.phi
  phi_linkobj <- make.link(phi_linkstr)
  phi_linkfun <- phi_linkobj$linkfun
  phi_linkinv <- phi_linkobj$linkinv
  phi_mu.eta <- phi_linkobj$mu.eta  

  ## y* transformation
  ystar <- qlogis(y)

  ## control parameters
  phi_full <- control$phi
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  ocontrol <- control
  control$phi <- control$method <- control$hessian <- control$start <- NULL

  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, linkfun(y), weights, offset = offset)
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    sigma2 <- sum(weights * res^2)/((sum(weights) - k) * (dlink)^2)
    phi_y <- weights * yhat * (1 - yhat)/(sum(weights) * sigma2) - 1/n
    phi <- rep(0, ncol(z))
    phi[1] <- phi_linkfun(sum(phi_y))
    ## i.e., start out from the fixed dispersion model as described
    ## in Ferrari & Cribari-Neto (2004) (and differing from Simas et al. 2009)
    ## An alternative would be
    ##   phi <- lm.wfit(z, phi_linkfun(phi_y), weights)$coefficients
    ## but that only works in general if all(phi_y > 0) which is not necessarily
    ## the case.
    start <- list(mean = beta, precision = phi)
  }
  if(is.list(start)) start <- do.call("c", start)

  ## objective function and gradient
  loglikfun <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    mu <- linkinv(x %*% beta + offset)
    phi <- phi_linkinv(z %*% gamma)
    if(any(!is.finite(phi))) NaN else { ## catch extreme cases without warning
      ## Use dbeta() instead of 'textbook' formula:
      ## ll <- lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
      ##   (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y)
      ll <- suppressWarnings(dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE))      
      if(any(!is.finite(ll))) NaN else sum(weights * ll) ## again: catch extreme cases without warning
    }
  }
  gradfun <- function(par) {
    beta <- par[1:k]
    gamma <- par[-(1:k)]
    eta <- as.vector(x %*% beta + offset)
    phi_eta <- as.vector(z %*% gamma)
    mu <- linkinv(eta)
    phi <- phi_linkinv(phi_eta)
    mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
    rval <- cbind(
      phi * (ystar - mustar) * mu.eta(eta) * weights * x,
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
        phi_mu.eta(phi_eta) * weights * z      
    )
    colSums(rval)
  }

  ## optimize likelihood  
  opt <- optim(par = start, fn = loglikfun, gr = gradfun,
    method = method, hessian = hessian, control = control)
  if(opt$convergence > 0) warning("optimization failed to converge")

  ## extract fitted values/parameters
  beta <- as.vector(opt$par[1:k])
  gamma <- as.vector(opt$par[-(1:k)])
  eta <- as.vector(x %*% beta + offset)
  phi_eta <- as.vector(z %*% gamma)
  mu <- linkinv(eta)
  phi <- phi_linkinv(phi_eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)
  pseudor2 <- if(var(eta) * var(ystar) <= 0) NA else cor(eta, linkfun(y))^2

  ## compute analytical covariance matrix
  ## <FIXME> 
  ## ## old code based on Ferrari & Cribari-Neto (2004)
  ## ## compute diagonal of T
  ## Tdiag <- mu.eta(eta)
  ## ## compute w
  ## w <- weights * phi * (psi1 + psi2) * Tdiag^2
  ## ## compute vector c
  ## vc <- phi * (psi1 * mu - psi2 * (1 - mu))
  ## ## compute d
  ## d <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
  ## ## compute (X'W X)^(-1)
  ## xwx1 <- chol2inv(qr.R(qr(sqrt(w) * x)))
  ## ## compute X'Tc
  ## xtc <- as.vector(t(x) %*% (Tdiag * vc))  
  ## ## compute gamma
  ## xwtc <- as.vector(xwx1 %*% xtc)
  ## gamma <- sum(d) - sum(xtc * xwtc)/phi  
  ## ## compute components of K^(-1)
  ## Kbb <- (xwx1/phi) %*% (diag(k) + outer(xtc, xwtc)/(gamma*phi))
  ## Kpp <- (1/gamma)
  ## Kbp <- -as.vector(xwx1 %*% xtc)/(gamma * phi)
  ## ## put together covariance matrix
  ## vcov <- rbind(cbind(Kbb, Kbp), c(Kbp, Kpp))
  ##
  ## instead: new code based on Simas et al. (2009): more general, possibly less exact
  ## due to solve(), maybe use inverse for partitioned matrices?
  ## </FIXME>
  ## auxiliary transformations
  a <- psi1 + psi2
  b <- psi1 * mu^2 + psi2 * (1-mu)^2 - trigamma(phi)
  ## compute elements of W
  wbb <- phi^2 * a * mu.eta(eta)^2
  wpp <- b * phi_mu.eta(phi_eta)^2
  wbp <- phi * (mu * a - psi2) * mu.eta(eta) * phi_mu.eta(phi_eta)
  ## compute elements of K
  kbb <- crossprod(sqrt(weights) * sqrt(wbb) * x)
  kpp <- crossprod(sqrt(weights) * sqrt(wpp) * z)
  kbp <- crossprod(weights * wbp * x, z)
  if(!hessian) {
    ## put together K
    hess <- cbind(rbind(kbb, t(kbp)), rbind(kbp, kpp))
    ## compute K^(-1)
    kbb1 <- chol2inv(qr.R(qr(sqrt(weights) * sqrt(wbb) * x)))
    kpp1 <- solve(kpp - t(kbp) %*% kbb1 %*% kbp)
    vcov <- cbind(rbind(kbb1 + kbb1 %*% kbp %*% kpp1 %*% t(kbp) %*% kbb1,
      -kpp1 %*% t(kbp) %*% kbb1), rbind(-kbb1 %*% kbp %*% kpp1, kpp1))
    ## vcov <- solve(hess)
  } else {
    ## if specified, use numerical Hessian instead of analytical solution
    hess <- -as.matrix(opt$hessian)
    vcov <- solve(hess)
  }

  ## names
  names(beta) <- colnames(x)
  names(gamma) <- if(phi_const & phi_linkstr == "identity") "(phi)" else colnames(z)
  rownames(vcov) <- colnames(vcov) <- rownames(hess) <- colnames(hess) <- c(colnames(x), 
    if(phi_const & phi_linkstr == "identity") "(phi)" else paste("(phi)", colnames(z), sep = "_"))

  ## set up return value
  rval <- list(  
    coefficients = list(mean = beta, precision = gamma),
    residuals = y - mu,
    fitted.values = structure(mu, .Names = names(y)),
    optim = opt,
    method = method,
    control = control,
    start = start,
    weights = if(identical(as.vector(weights), rep(1, n))) NULL else weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    df.null = n - 2,
    df.residual = n - k - m,
    phi = phi_full,
    loglik = opt$value,
    vcov = vcov,
    hessian = hess,
    pseudo.r.squared = pseudor2,
    link = list(mean = linkobj, precision = phi_linkobj),
    converged = opt$convergence < 1    
  )

  return(rval)
}

print.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if(length(coef(x))) {
      cat(paste("Coefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
      print.default(format(x$coefficients$mean, digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
      if(x$phi) {
        cat(paste("Phi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
        print.default(format(x$coefficients$precision, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      }
    }
    else cat("No coefficients\n\n")
  }
  
  invisible(x)
}

summary.betareg <- function(object, phi = NULL, type = "sweighted2", ...)
{
  ## treat phi as full model parameter?
  if(!is.null(phi)) object$phi <- phi
  
  ## residuals
  type <- match.arg(type, c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2"))
  object$residuals <- residuals(object, type = type)
  object$residuals.type <- type
  
  ## extend coefficient table
  k <- length(object$coefficients$mean)
  cf <- as.vector(do.call("c", object$coefficients))
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(mean = cf[1:k, , drop = FALSE], precision = cf[-(1:k), , drop = FALSE])
  rownames(cf$mean) <- names(object$coefficients$mean)
  rownames(cf$precision) <- names(object$coefficients$precision)
  object$coefficients <- cf
  
  ## number of iterations
  object$iterations <- as.vector(tail(na.omit(object$optim$count), 1))  
  
  ## delete some slots
  object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- NULL

  ## return
  class(object) <- "summary.betareg"
  object
}

print.summary.betareg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    types <- c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
    Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals",
      "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")  
    cat(sprintf("%s:\n", Types[types == match.arg(x$residuals.type, types)]))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  
    cat(paste("\nCoefficients (mean model with ", x$link$mean$name, " link):\n", sep = ""))
    printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
  
    if(x$phi) {
      cat(paste("\nPhi coefficients (precision model with ", x$link$precision$name, " link):\n", sep = ""))
      printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
    }
    
    if(getOption("show.signif.stars") & any(do.call("rbind", x$coefficients)[,4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
  
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, NROW)), "Df")
    if(!is.na(x$pseudo.r.squared)) cat("\nPseudo R-squared:", formatC(x$pseudo.r.squared, digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
  }
  
  invisible(x)
}

predict.betareg <- function(object, newdata = NULL,
  type = c("response", "link", "precision", "variance"), na.action = na.pass, ...) 
{
  type <- match.arg(type)
  
  if(missing(newdata)) {

    rval <- switch(type,
      "response" = {
        object$fitted.values
      },
      "link" = {
        object$link$mean$linkfun(object$fitted.values)
      },      
      "precision" = {
        gamma <- object$coefficients$precision
        z <- if(is.null(object$x)) model.matrix(object, model = "precision") else object$x$precision
	object$link$precision$linkinv(drop(z %*% gamma))
      },
      "variance" = {
        gamma <- object$coefficients$precision
        z <- if(is.null(object$x)) model.matrix(object, model = "precision") else object$x$precision
	phi <- object$link$precision$linkinv(drop(z %*% gamma))
	object$fitted.values * (1 - object$fitted.values) / (1 + phi)
      }
    )
    return(rval)

  } else {

    tnam <- switch(type, 
      "response" = "mean",
      "link" = "mean",
      "precision" = "precision",
      "variance" = "full")
      
    mf <- model.frame(delete.response(object$terms[[tnam]]), newdata, na.action = na.action, xlev = object$levels[[tnam]])
    if(type %in% c("response", "link", "variance")) {
      X <- model.matrix(delete.response(object$terms$mean), mf, contrasts = object$contrasts$mean)
      offset <- if(!is.null(off.num <- attr(object$terms$full, "offset")))
    	eval(attr(object$terms$full, "variables")[[off.num + 1]], newdata)
        else if(!is.null(object$offset)) eval(object$call$offset, newdata)
      if(is.null(offset)) offset <- rep(0, NROW(X))
    }
    if(type %in% c("precision", "variance")) {
      Z <- model.matrix(object$terms$precision, mf, contrasts = object$contrasts$precision)
    }

    rval <- switch(type,    
      "response" = {
        object$link$mean$linkinv(drop(X %*% object$coefficients$mean + offset))
      },      
      "link" = {
        drop(X %*% object$coefficients$mean + offset)
      },      
      "precision" = {
        object$link$precision$linkinv(drop(Z %*% object$coefficients$precision))
      },
      "variance" = {
        mu <- object$link$mean$linkinv(drop(X %*% object$coefficients$mean + offset))
        phi <- object$link$precision$linkinv(drop(Z %*% object$coefficients$precision))
	mu * (1 - mu) / (1 + phi)
      }
    )
    return(rval)

  }
}

coef.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...) {
  cf <- object$coefficients

  model <- if(is.null(phi)) {
    if(missing(model)) ifelse(object$phi, "full", "mean") else match.arg(model)
  } else {
    if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
    ifelse(phi, "full", "mean")
  }
  
  switch(model,
    "mean" = {
      cf$mean
    },
    "precision" = {
      cf$precision
    },
    "full" = {
      nam1 <- names(cf$mean)
      nam2 <- names(cf$precision)
      cf <- c(cf$mean, cf$precision)
      names(cf) <- c(nam1, if(identical(nam2, "(phi)")) "(phi)" else paste("(phi)", nam2, sep = "_"))
      cf
    }
  )
}

vcov.betareg <- function(object, model = c("full", "mean", "precision"), phi = NULL, ...) {
  vc <- object$vcov
  k <- length(object$coefficients$mean)

  model <- if(is.null(phi)) {
    if(missing(model)) ifelse(object$phi, "full", "mean") else match.arg(model)
  } else {
    if(!missing(model)) warning("only one of 'model' and 'phi' should be specified: 'model' ignored")
    ifelse(phi, "full", "mean")
  }
  
  switch(model,
    "mean" = {
      vc[1:k, 1:k, drop = FALSE]
    },
    "precision" = {
      vc <- vc[-(1:k), -(1:k), drop = FALSE]
      colnames(vc) <- rownames(vc) <- names(object$coefficients$precision)
      vc
    },
    "full" = {
      vc
    }
  )
}

bread.betareg <- function(x, phi = NULL, ...) {
  vcov(x, phi = phi) * x$n
}

estfun.betareg <- function(x, phi = NULL, ...)
{
  ## extract response y and regressors X and Z
  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  xmat <- if(is.null(x$x)) model.matrix(x, model = "mean") else x$x$mean
  zmat <- if(is.null(x$x)) model.matrix(x, model = "precision") else x$x$precision
  offset <- if(is.null(x$offset)) rep(0, NROW(xmat)) else x$offset
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  phi_full <- if(is.null(phi)) x$phi else phi
  
  ## extract coefficients
  beta <- x$coefficients$mean
  gamma <- x$coefficients$precision

  ## compute y*
  ystar <- qlogis(y)
  
  ## compute mu*
  eta <- as.vector(xmat %*% beta + offset)
  phi_eta <- as.vector(zmat %*% gamma)
  mu <- x$link$mean$linkinv(eta)
  phi <- x$link$precision$linkinv(phi_eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ## compute scores of beta
  rval <- phi * (ystar - mustar) * as.vector(x$link$mean$mu.eta(eta)) * wts * xmat

  ## combine with scores of phi
  if(phi_full) {
    rval <- cbind(rval,
      (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)) *
      as.vector(x$link$precision$mu.eta(phi_eta)) * wts * zmat)
    colnames(rval) <- names(coef(x, phi = phi_full))
  }

  attr(rval, "assign") <- NULL
  return(rval)
}

coeftest.betareg <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)  

logLik.betareg <- function(object, ...) {
  structure(object$loglik, df = sum(sapply(object$coefficients, length)), class = "logLik")
}

terms.betareg <- function(x, model = c("mean", "precision"), ...) {
  x$terms[[match.arg(model)]]
}

model.frame.betareg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  if(is.Formula(formula$formula)) formula$call$formula <- formula$formula <-
    formula(formula$formula, collapse = TRUE)
  formula$terms <- formula$terms$full
  NextMethod()
}

model.matrix.betareg <- function(object, model = c("mean", "precision"), ...) {
  model <- match.arg(model)
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
    else model.matrix(object$terms[[model]], model.frame(object), contrasts = object$contrasts[[model]])
  return(rval)
}

residuals.betareg <- function(object,
  type = c("sweighted2", "deviance", "pearson", "response", "weighted", "sweighted"), ...)
{
  ## raw response residuals and desired type
  res <- object$residuals
  type <- match.arg(type)
  if(type == "response") return(res)

  ## extract fitted information
  y <- if(is.null(object$y)) model.response(model.frame(object)) else object$y
  mu <- fitted(object)
  wts <- weights(object)
  if(is.null(wts)) wts <- 1
  phi <- predict(object, type = "precision")
  
  res <- switch(type,
  
    "pearson" = {
      sqrt(wts) * res / sqrt(mu * (1 - mu) / (1 + phi))
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
    h <- hatvalues(model)
    k <- length(model$coefficients$mean)
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
