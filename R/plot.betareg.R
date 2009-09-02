plot.betareg <- function(x, which = 1:4,
  caption = c("Residuals vs indices of obs.", "Cook's distance plot",
    "Generalized leverage vs predicted values", "Residuals vs linear predictor", 
    "Half-normal plot of residuals"),
    sub.caption = paste(deparse(x$call), collapse = "\n"), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., type = "deviance", nsim = 100, level = 0.9)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("`which' must be in 1:5")
    
  types <- c("pearson", "deviance", "response", "weighted", "sweighted", "sweighted2")
  Types <- c("Pearson residuals", "Deviance residuals", "Raw response residuals",
    "Weighted residuals", "Standardized weighted residuals", "Standardized weighted residuals 2")
  type <- match.arg(type, types)
  Type <- Types[type == types]

  res <- residuals(x, type = type)
  n <- length(res)
  k <- length(x$coefficients) - 1
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if(show[1]) {
    plot(1:n, res, xlab = "Obs. number", ylab = Type, main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[2]) {
    plot(1:n, cooks.distance(x),
      xlab = "Obs. number", ylab = "Cook's distance", type = "h", main = main)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
  }
  if(show[3]) {
    plot(fitted(x), gleverage(x),
      xlab = "Predicted values", ylab = "Generalized leverage", main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25)
  }
  if(show[4]) {
    plot(predict(x, type = "link"), res,
      xlab = "Linear predictor", ylab = Type, main = main, ...)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[4], 3, 0.25)
    abline(h = 0, lty = 3, col = "gray")
  }
  if(show[5]) {
    hn <- halfnormal.betareg(x, nsim = nsim, level = level, type = type)
    plot(hn[,1], hn[,2], ylim = range(hn[,-1]), main = main,
      xlab = "Normal quantiles", ylab = paste(Type, "(absolute values)"), ...)
    lines(hn[,1], hn[,3],lty = 2)
    lines(hn[,1], hn[,4],lty = 1)
    lines(hn[,1], hn[,5],lty = 1)
    if(one.fig) title(sub = sub.caption, ...)
    mtext(caption[5], 3, 0.25)
  }

  if(!one.fig && par("oma")[3] >= 1) mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}

halfnormal.betareg <- function(model, nsim = 100, level = 0.90, type = "deviance")
{
  ## extract response y and regressors X
  y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
  x <- if(is.null(model$x)) model.matrix(model) else model$x
  offset <- if(is.null(model$offset)) rep(0, NROW(x)) else model$offset
  wts <- weights(model)

  n <- NROW(x)
  alpha <- (1 - level)/2
  mu <- fitted(model)
  phi <- tail(model$coefficients, 1)    
  res <- residuals(model, type = type)

  e <- matrix(0, n, nsim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:nsim) {
    ysim <- rbeta(n, mu * phi, (1 - mu) * phi)
    fit <- betareg(ysim ~ 0 + x, weights = wts, offset = offset)
    e[,i] <- sort(abs(residuals(fit, type = type)))
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo, alpha)
    e2[i] <- quantile(eo, 1 - alpha)
  }
  
  e0 <- apply(e, 1, median)
  qq <- qnorm((n + 1:n + 0.5)/(2 * n + 1.125))
  
  cbind(qq, sort(abs(res)), e0, e1, e2)  
}
