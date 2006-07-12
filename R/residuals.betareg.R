"residuals.betareg" <-
function(object, type = c("standardized","raw","deviance"),...)
{
    type <- match.arg(type)
    muhat <- object$fitted
    r = object$y - muhat
    rd <- object$resd
    rd = as.vector(rd)
    phihat <- c(object$coef[object$k+1])
    sqrtvar <- sqrt(muhat*(1-muhat)/(1+phihat))
    r = as.vector(r)
    res <- switch(type, standardized = as.vector(r/sqrtvar),usual=r,deviance = rd)
        if (is.null(object$na.action)) 
        res
    else naresid(object$na.action, as.vector(r/sqrtvar))
}
