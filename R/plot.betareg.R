"plot.betareg" <-
function (x, which = 1:4, caption = c("Deviance residuals vs indices of obs.", 
    "Standardized residuals vs indices of obs.", "Generalized leverage vs. Predicted values", "Cook's distance plot"), 
    panel = points, sub.caption = deparse(x$call), main = "", 
    ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75) 
{
    if (!inherits(x, "betareg")) 
        stop("Use only with 'betareg' objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
    rdev <- residuals(x,type="deviance")
    n <- length(rdev)
    h <- x$h
    k <- x$k
    gl <- x$GL
    gl <- diag(gl)
    yh <- predict(x)
    sr <- residuals(x)
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
        plot(rdev, xlab = "Indices of obs.", ylab = "Deviance residuals", main = main, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[1], 3, 0.25)
        abline(h = 0, lty = 3, col = "gray")
    }
    if (show[2]) {
        plot(sr, xlab = "Obs. number", ylab = "Standardized residuals", main = main, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[2], 3, 0.25)
        abline(h = 0, lty = 3, col = "gray")
    }
    if (show[3]) {
	        plot(yh, gl, xlab = "Fitted values", ylab = "Generalized leverage", main = main,...)
        panel(yh, gl, ...)
        if (one.fig) 
            title(sub = sub.caption, ...)
        mtext(caption[3], 3, 0.25)
    }
    if (show[4]) {
	plot(h*(sr^2)/(k*(1-h)^2),xlab="Obs. number",ylab = "Cook's distance",type="h",main=main)
	mtext(caption[4], 3, 0.25)
       if (one.fig) 
            title(sub = sub.caption, ...)
    }
    if (!one.fig && par("oma")[3] >= 1) 
        mtext(sub.caption, outer = TRUE, cex = 1.25)
    invisible()
}
