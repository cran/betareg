"predict.betareg" <-
function(object, terms = object$x[,-1],... ) 
{
   coef <- (object$coeff)[1:object$k]
   x <- terms 
   x <- matrix(x,length(x)/(ncol(object$x)-1),ncol(object$x)-1)
   x <- cbind(1,x)
   if(ncol(x) != ncol(object$x))
       stop("The number of columns must be equal to the number of coefficients.")
   pred = x%*%coef
   pred = object$linkinv(pred)
   pred = as.vector(pred)
   pred
}

