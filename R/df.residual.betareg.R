df.residual.betareg <- function (object, ...)
{
   df <- length(object$y) - object$k
   df
}
