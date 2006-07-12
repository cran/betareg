"summary.betareg" <-
function (object,...) 
{
    z <- object
    ans <<- z[c("call", "terms")]
    class(ans) <<- "summary.betareg"
    y <- object$y
    nulldev <- z$nulldev
    resdev <- sum( (residuals(z,type="deviance"))^2)
    ans$coefficients <<- z$coef
    ans$std <<- z$stder
    ans$zstats <<- z$zstats
    ans$pvalues <<- z$pvalues
    k <- z$k
    y <- z$y
    lik1 <- logLik(z)
    lik2 <- logLik(betareg(y ~ 1))
    likratio <- 2*(lik1 - lik2)
    ans$est <<- cbind(format(ans$coeff[1:k],digits=5),format(ans$std[1:k],digits=4,sci = TRUE),format(ans$zstats[1:k],digits=3),format(ans$pvalues[1:k],digits=3))
    dimnames(ans$est) <- list(NULL, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    rownames(ans$est) <- names(ans$coef[-(k+1)])
    ans$Pseudo.R2 <- z$Pseudo.R2
    cat("\nCall:\n", deparse(ans$call), "\n\n", sep = "")
    if (length(coef(ans))) {
        cat("Deviance Residuals:\n")
        print(summary(residuals(z,type="deviance"),digits=6)[-4])
        cat("\n")
        cat("Coefficients:\n")
        print.default(format(ans$est),print.gap = 2, quote = FALSE)
        cat("---\n")
        cat("\n")
        cat("Estimated precision parameter (phi): ")
        cat(ans$coefficients[k+1])
        cat(" with s.e. ")
        cat(ans$std[k+1])
        cat("\n\n")
        cat("    Null Deviance: ")
        cat(nulldev)
        cat(" on ")
        cat(length(y)-1)
        cat(" degrees of freedom\n")
        cat("Residual Deviance: ")
        cat(resdev)
        cat(" on ")
        cat(length(y)-k)
        cat(" degrees of freedom\n")
        cat("Log-Likelihood Ratio Statistic: ")
        cat(likratio)
        cat(" on ")
        cat(k-1)
        cat(" degrees of freedom\n\n")
        cat("Pseudo R^2: ")
        cat(format(ans$Pseudo.R2))
        cat("\t\t AIC: ")
        cat(AIC(z))
        cat("\n\n")
    }
    invisible(ans)
}
