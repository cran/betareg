gen.lev.betareg <- function(x, ...)
{
    if (!inherits(x, "betareg")) 
    stop("Use only with 'betareg' objects")
    gl <- x$GL
    gl <- diag(gl)
    gl
}