anova.betareg <- function(object, object2, ...) {
if(missing(object2)){
cat("Analysis of Variance Table\n\n")
cat("Model: ")
print(object$formula)
y2 <- object$y
obj2 <- betareg(y2 ~ 1)
df <- c(object$k,obj2$k)
loglik <- c(logLik(object),logLik(obj2))
diflik <- logLik(object)-logLik(obj2)
difdf <- object$k - obj2$k
diflik <- 2*diflik
pval <- 1-pchisq(diflik,difdf)
est <<- cbind(format(df,digits=5),format(loglik,digits=6),c("",format(diflik,digits=5)),c("",format(pval,digits=5)))
dimnames(est) <- list(NULL, c("df", "Log. Lik", "Ratio", "P(>|Chi|)"))
rownames(est) <- c("1","2")
print.default(format(est),print.gap = 2, quote = FALSE)
}
else{
cat("Analysis of Variance Table\n\n")
cat("Model 1: ")
print(object$formula)
cat("Model 2: ")
print(object2$formula)
df <- c(object$k,object2$k)
loglik <- c(logLik(object),logLik(object2))
diflik <- -1* (logLik(object)-logLik(object2))
difdf <- -1 * (object$k - object2$k)
diflik <- 2*diflik
pval <- 1-pchisq(diflik,difdf)
est <<- cbind(format(df,digits=5),format(loglik,digits=6),c("",format(diflik,digits=5)),c("",format(pval,digits=5)))
dimnames(est) <- list(NULL, c("df", "Log. Lik", "Ratio", "P(>|Chi|)"))
rownames(est) <- c("1","2")
print.default(format(est),print.gap = 2, quote = FALSE)
}
}