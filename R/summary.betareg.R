"summary.betareg" <-
function (object,...) 
{
    z <- object
    ans <<- z[c("call", "terms")]
    class(ans) <<- "summary.betareg"
    ans$coefficients <<- z$coef
    ans$std <<- z$stder
    ans$zstats <<- z$zstats
    ans$pvalues <<- z$pvalues
    k <- z$k
    ans$est <<- cbind(format(ans$coeff),format(ans$std),c(format(ans$zstats[1:k]),""),c(format(ans$pvalues[1:k]),""))
    dimnames(ans$est) <- list(NULL, c("estimates", "std. errors","z-stats", "p-value"))
    rownames(ans$est) <- names(ans$coef)
    ans$Pseudo.R2 <- z$Pseudo.R2
    cat("\nCall:\n", deparse(ans$call), "\n\n", sep = "")
    if (length(coef(ans))) {
        cat("Coefficients:\n")
        print.default(format(ans$est),print.gap = 2, quote = FALSE)
        cat("\n")
        cat("Pseudo R^2:")
        print.default(format(ans$Pseudo.R2), print.gap = 2, quote = FALSE)
    }
    invisible(ans)
}
