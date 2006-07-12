logLik.betareg <- function(object, ...)
{
    val <- object$value
    p <- object$k
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}
