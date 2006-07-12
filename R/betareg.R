"betareg" <- function (formula, link = "logit", data)
{
  call <- match.call()
  m <- match.call(expand = FALSE)
  link = c(link)
  m$link <- NULL
  m$x <- m$y
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.frame(sys.parent()))
  Terms <- attr(m, "terms")
  Y <- model.extract(m, response)
  X <- model.matrix(Terms, m, contrasts)
  if (min(Y) <= 0 || max(Y) >= 1)
    stop("OUT OF RANGE (0,1)!")
  offset <- model.offset(m)
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  link1 <- structure(list(link = linktemp, linkfun = stats$linkfun,
                          linkinv = stats$linkinv, mu.eta = stats$mu.eta, diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t))),
                          T = function(etahat) diag(c(stats$mu.eta(etahat)))))
  fit <- c()
  fit1 <- br.fit(X, Y, link1)
  fit$call <- call
  fit$funlink <- link1
  fit$linkinv <- link1$linkinv
  fit$coefficients <- c(fit1$coeff, structure(fit1$phiest,
                                              .Names = c("phi")))
  fit$stder <- c(fit1$stder, fit1$phist)
  fit$zstats <- fit1$zstats
  fit$pvalues <- fit1$pvalues
  fit$fitted.values <- fit1$fitted
  fit$linear.predictor <- fit1$linpred
  fit$residuals <- fit1$res
  fit$k <- fit1$k
  fit$nulldev <- fit1$nulldev
  fit$value <- fit1$value
  fit$h <- fit1$h
  fit$GL <- fit1$GL
  fit$terms <- Terms
  fit$x <- X
  fit$y <- Y
  fit$resd <- fit1$resd
  fit$df.residual <- fit1$df.residual
  fit$resstd <- fit1$resstd
  fit$Pseudo.R2 <- cbind(fit1$pseudor2)
     fit$etahat <- fit1$etahat
     fit$sigma2 <- fit1$sigma2
     fit$phi <- fit1$phi
  fit$formula <- formula
  colnames(fit$Pseudo.R2) <- c("")
  rownames(fit$Pseudo.R2) <- c("")
  attr(fit, "na.message") <- attr(m, "na.message")
  class(fit) <- c("betareg", "lm")
  fit
}
