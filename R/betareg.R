betareg <- function(formula, data, subset, na.action, weights, offset,
                    link = c("logit", "probit", "cloglog"),
 		    control = betareg.control(...),
		    model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## evaluate model.frame
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrix, response
  mt <- terms(formula, data = data)
  X <- model.matrix(mt, mf)
  Y <- model.response(mf, "numeric")

  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(!(min(Y) >= 0 & max(Y) <= 1)) stop("invalid dependent variable, must be in [0, 1]")

  ## convenience variables
  n <- length(Y)
  k <- NCOL(X)

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
  rval <- betareg.fit(X, Y, weights, offset, match.arg(link), control)

  ## further model information
  rval$call <- cl
  rval$formula <- formula
  rval$terms <- mt
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- attr(X, "contrasts")
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- X

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

betareg.fit <- function(x, y, weights = NULL, offset = NULL,
  link = "logit", control = betareg.control())
{
  ## link processing
  linkstr <- link
  linkobj <- make.link(linkstr)
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta

  ## response and regressor matrix
  n <- NROW(x)
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep(1, n)
  if(is.null(offset)) offset <- rep(0, n)
  ystar <- linkfun(y)

  ## control parameters
  phi_full <- control$phi
  method <- control$method
  hessian <- control$hessian
  start <- control$start
  ocontrol <- control
  control$phi <- control$method <- control$hessian <- control$start <- NULL

  ## starting values
  if(is.null(start)) {
    auxreg <- lm.wfit(x, ystar, weights, offset = offset)
    beta <- auxreg$coefficients
    yhat <- linkinv(auxreg$fitted.values)
    dlink <- 1/mu.eta(linkfun(yhat))
    res <- auxreg$residuals
    sigma2 <- sum((weights * res)^2)/((n - k) * (dlink)^2)
    phi <- mean(yhat * (1 - yhat)/sigma2 - 1)
    start <- c(beta, phi)
  }

  ## objective function and gradient
  loglikfun <- function(par) {
    beta <- par[1:k]
    phi <- par[k+1]
    mu <- linkinv(x %*% beta + offset)
    ll <- lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) + 
      (mu * phi - 1) * log(y) + ((1 - mu) * phi - 1) * log(1 - y)
    sum(weights * ll)
  }
  gradfun <- function(par) {
    beta <- par[1:k]
    phi <- par[k+1]
    eta <- as.vector(x %*% beta + offset)
    mu <- linkinv(eta)
    mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
    rval <- cbind(phi * (ystar - mustar) * mu.eta(eta) * weights * x,
      weights * (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)))
    colSums(rval)
  }

  ## optimize likelihood  
  opt <- optim(par = start, fn = loglikfun, gr = gradfun,
    method = method, hessian = hessian, control = control)
  if(opt$convergence > 0) warning("optimization failed to converge")

  ## extract fitted values/parameters
  beta <- as.vector(opt$par[1:k])
  phi <- as.vector(opt$par[k+1])
  eta <- as.vector(x %*% beta + offset)
  mu <- linkinv(eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
  psi1 <- trigamma(mu * phi)
  psi2 <- trigamma((1 - mu) * phi)
  pseudor2 <- if(var(eta) * var(ystar) <= 0) NA else cor(eta, ystar)^2

  ## compute analytical covariance matrix
  ## compute diagonal of T
  Tdiag <- mu.eta(eta)
  ## compute w
  w <- weights * phi * (psi1 + psi2) * Tdiag^2
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
  ## compute components of K^(-1)
  Kbb <- (xwx1/phi) %*% (diag(k) + outer(xtc, xwtc)/(gamma*phi))
  Kpp <- (1/gamma)
  Kbp <- -as.vector(xwx1 %*% xtc)/(gamma * phi)
  ## put together covariance matrix
  vcov <- rbind(cbind(Kbb, Kbp), c(Kbp, Kpp))
  rownames(vcov) <- colnames(vcov) <- c(colnames(x), "(phi)")

  ## if specified, use Hessian instead of analytical solution
  if(hessian) vcov <- -solve(as.matrix(opt$hessian))

  ## set up return value
  rval <- list(  
    coefficients = structure(as.vector(opt$par), .Names = c(colnames(x), "(phi)")),
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
    df.residual = n - k,
    phi = phi_full,
    loglik = opt$value,
    vcov = vcov,
    pseudo.R.squared = pseudor2,
    link = linkobj,
    converged = opt$convergence < 1    
  )

  return(rval)
}

print.betareg <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
        cat("\n")
	if(!x$phi) cat(sprintf("Estimated precision parameter (phi): %s\n\n",
	  format(tail(x$coefficients, 1), digits = digits)))
    }
    else cat("No coefficients\n\n")
  }
  
  invisible(x)
}

summary.betareg <- function(object, phi = NULL, ...)
{
  ## treat phi as full model parameter?
  if(!is.null(phi)) object$phi <- phi
  
  ## store deviance residuals
  object$residuals <- residuals(object, type = "deviance")
  
  ## extend coefficient table
  cf <- object$coefficients
  se <- sqrt(diag(object$vcov))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf
  
  ## number of iterations
  object$iterations <- tail(na.omit(object$optim$count), 1)
  
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
    cat("Deviance residuals:\n")
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
      .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  
    cf <- x$coefficients
    if(!x$phi) {
      phi <- cf[NROW(cf),]
      cf <- cf[-NROW(cf),, drop = FALSE]
    }
    cat("\nCoefficients:\n")
    printCoefmat(cf, digits = digits, ...)
    cat("\n")
  
    if(!x$phi) cat(sprintf("Estimated precision parameter (phi): %s\n",
      format(phi[1], digits = digits)))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual + 1, "Df\n")
    cat("Pseudo R-squared:", formatC(x$pseudo.R.squared, digits = digits))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
  }
  
  invisible(x)
}

predict.betareg <- function(object, newdata = NULL,
  type = c("response", "link"), na.action = na.pass, ...) 
{
  type <- match.arg(type)
  
  if(missing(newdata)) {
    if(type == "response") return(object$fitted.values)
      else return(object$link$linkfun(object$fitted.values))
  } else {
    mf <- model.frame(delete.response(object$terms), newdata, na.action = na.action, xlev = object$levels)
    X <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
    offset <- if(!is.null(off.num <- attr(object$terms, "offset")))
  	eval(attr(object$terms, "variables")[[off.num + 1]], newdata)
      else if(!is.null(object$offset)) eval(object$call$offset, newdata)
    if(is.null(offset)) offset <- rep(0, NROW(X))

    pred <- drop(X %*% coef(object, phi = FALSE) + offset)
    if(type == "response") return(object$link$linkinv(pred))
      else return(pred)
  }
}

coef.betareg <- function(object, phi = NULL, ...) {
  cf <- object$coefficients
  if(is.null(phi)) phi <- object$phi
  if(phi) cf else cf[-length(cf)]
}

vcov.betareg <- function(object, phi = NULL, ...) {
  vc <- object$vcov
  if(is.null(phi)) phi <- object$phi
  if(phi) vc else vc[-NROW(vc), -NROW(vc), drop = FALSE]
}

bread.betareg <- function(x, phi = NULL, ...) {
  vcov(x, phi = phi) * x$n
}

estfun.betareg <- function(x, phi = NULL, ...)
{
  ## extract response y and regressors X
  y <- if(is.null(x$y)) model.response(model.frame(x)) else x$y
  xmat <- if(is.null(x$x)) model.matrix(x) else x$x
  offset <- if(is.null(x$offset)) rep(0, NROW(xmat)) else x$offset
  wts <- weights(x)
  if(is.null(wts)) wts <- 1
  phi_full <- if(is.null(phi)) x$phi else phi
  
  ## extract coefficients
  beta <- x$coefficients
  phi <- beta[length(beta)]
  beta <- beta[-length(beta)]

  ## compute y*
  ystar <- as.vector(x$link$linkfun(y))
  
  ## compute mu*
  eta <- as.vector(xmat %*% beta + offset)
  mu <- x$link$linkinv(eta)
  mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)

  ## compute diagonal of matrix T
  Tdiag <- as.vector(x$link$mu.eta(eta))

  ## compute scores of beta
  rval <- phi * (ystar - mustar) * Tdiag * wts * xmat

  ## combine with scores of phi
  if(phi_full) rval <- cbind(rval,
    "(phi)" = wts * (mu * (ystar - mustar) + log(1-y) - digamma((1-mu)*phi) + digamma(phi)))

  attr(rval, "assign") <- NULL
  return(rval)
}

coeftest.betareg <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)  

logLik.betareg <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual + 1, class = "logLik")
}

terms.betareg <- function(x, ...) {
  x$terms
}

model.frame.betareg <- function(formula, ...) {
  if(!is.null(formula$model)) return(formula$model)
  NextMethod()
}

model.matrix.betareg <- function(object, ...) {
  rval <- if(!is.null(object$x)) object$x
    else model.matrix(object$terms, model.frame(object), contrasts = object$contrasts)
  return(rval)
}

residuals.betareg <- function(object,
  type = c("deviance", "pearson", "response", "weighted", "sweighted", "sweighted2"), ...)
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
  phi <- tail(object$coefficients, 1)
  
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
      ystar <- object$link$linkfun(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(phi * v)
    },
    
    "sweighted" = {
      ystar <- object$link$linkfun(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(v)
    },
    
    "sweighted2" = {
      ystar <- object$link$linkfun(y)
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      v <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      sqrt(wts) * (ystar - mustar) / sqrt(v * hatvalues(object))
    })

  return(res)
}

cooks.distance.betareg <- function (model, ...) 
{
    h <- hatvalues(model)
    k <- length(model$coefficients) - 1
    res <- residuals(model, type = "pearson")
    h * (res^2)/(k * (1 - h)^2)
}
