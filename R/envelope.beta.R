envelope.beta <- function(model=fit.model,sim=100,conf=.90, pch="+",font.main=1, cex.main=1.5, type = c("standardized","deviance")) {
   type <- match.arg(type)
  if(!any(class(model)=="betareg"))
  {stop("The model must be from the class betareg")}
  main = switch(type, standardized = "Half-Normal Plot of Standardized Residuals", deviance = "Half-Normal Plot of Deviance Residuals")
  ylab= switch(type, standardized = "Absolute Values of Standardized Residuals", deviance = "Absolute Values of Deviance Residuals")
  xlab="Half-Normal Quantile"
  alfa <-(1-conf)/2
  X <- model$x
  y <-model$y
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  m <- model$fitted
  
  
###The line below is to avoid division by 0 when studentize the residuals, but trying to keep the leverage value high.
  h[round(h,15)==1]<-0.999999999999999
  
  si <- model$sigma2
  r <- model$residuals
  res <- switch(type, standardized = model$resstd, deviance = model$resd)
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  phi <- model$phi
  
  for(i in 1:sim) {
    resp <- rbeta(n, m*phi, (1-m)*phi)
    fit <- betareg(resp~X-1)
    ti <- fit$residuals/(model$sigma2*sqrt(1-h))
    eo <- switch(type,standardized = sort(abs(fit$resstd)), deviance = sort(abs(fit$resd)))
    e[,i] <- eo
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa)
    e2[i] <- quantile(eo,1-alfa)
  }
  
  med <- apply(e,1,median)
  qq <- qnorm((n+1:n+.5)/(2*n+1.125))
  plot(qq, sort(abs(res)), ylim=range(abs(res), e1, e2), pch=pch,
       main=main, xlab=xlab, ylab=ylab, cex.main=cex.main, font.main=font.main)
  
  lines(qq,e1,lty=1)
  lines(qq,e2,lty=1)
  lines(qq,med,lty=2) 
}

																												
