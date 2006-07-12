"predict.betareg" <-
function(object, newdata = NULL, type = c("link", "response"), ... ) 
{
   type <- match.arg(type)
         if (missing(newdata)) {
             pred <- switch(type, link = object$linear.predictor, response = object$fitted.values)
                             }
         else {
            dados <- model.frame(newdata)
            dados <- as.matrix(dados)
            coef <- (object$coeff)[1:object$k]
            x <- cbind(1,dados)
            pred = x%*%coef
            switch(type, response = {
                pred <- object$linkinv(pred)
            }, link = , terms = )
              }
pred <- as.vector(pred)
pred
}
