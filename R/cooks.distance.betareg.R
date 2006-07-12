cooks.distance.betareg <- function(model, ...)
{
    if (!inherits(model, "betareg")) 
        stop("Use only with 'betareg' objects")
    h <- model$h
    k <- model$k
    sr <- residuals(model)
    cook <- h*(sr^2)/(k*(1-h)^2)
    cook
}