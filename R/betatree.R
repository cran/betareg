## high-level convenience interface
betatree <- function(formula, partition, data, subset = NULL, na.action = na.omit, weights, offset, cluster,
		    link = "logit", link.phi = "log", control = betareg.control(), ...)
{
  ## use dots for setting up mob_control
  control <- partykit::mob_control(...)
  control$xtype <- control$ytype <- "data.frame"

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula to three parts if necessary
  f <- if(missing(partition)) Formula::Formula(formula) else Formula::as.Formula(formula, partition)
  if(length(f)[2L] == 1L) {
    attr(f, "rhs") <- c(list(1), list(1), attr(f, "rhs"))
    formula <- Formula::as.Formula(formula(f))
  } else if(length(f)[2L] == 2L) {
    attr(f, "rhs") <- c(attr(f, "rhs")[[1L]], list(1), attr(f, "rhs")[[2L]])
    formula <- Formula::as.Formula(formula(f))
  } else {
    formula <- f
  }
  ## formula for mob needs to collapse first two parts
  mob_formula <- formula(Formula::as.Formula(
    formula(formula, rhs = 1L:2L, collapse = TRUE),
    formula(formula, lhs = 0L, rhs = 3L)
  ))
  ## call for beta regression
  br_call <- match.call(expand.dots = FALSE)
  br_call$partition <- br_call$cluster <- br_call[["..."]] <- NULL
  br_call$formula <- formula(formula, lhs = 1L, rhs = 1L:2L)
  
  ## terms
  ft <- terms(formula, data = data, lhs = 1L, rhs = 1L:2L)
  xt <- terms(formula, data = data, lhs = 1L, rhs = 1L)
  zt <- terms(formula, data = data, lhs = 0L, rhs = 2L)
  
  ## actual fitting function for mob()
  betafit <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
    estfun = FALSE, object = FALSE)
  {
    ## catch control arguments
    args <- list(...)
    ctrl <- list(start = start)
    anam <- names(args)
    anam <- anam[!(anam %in% c("link", "link.phi", "type"))]
    for(n in anam) {
      ctrl[[n]] <- args[[n]]
      args[[n]] <- NULL
    }
    args$control <- do.call("betareg.control", ctrl)
  
    ## extract response and regressors
    mf <- cbind(y, x)
    attr(mf, "terms") <- ft
    y <- y[[1L]]
    xx <- model.matrix(xt, mf)
    xz <- model.matrix(zt, mf)    
  
    ## call betareg fitting function
    args <- c(list(x = xx, y = y, z = xz, weights = weights, offset = offset), args)
    obj <- do.call("betareg.fit", args)

    ## list structure
    rval <- list(
      coefficients = coef.betareg(obj),
      objfun = obj$loglik,
      estfun = NULL,
      object = NULL
    )

    ## add model (if desired)
    if(estfun | object) {
      class(obj) <- "betareg"
      obj$contrasts <- attr(x, "contrasts")
      obj$xlevels <- attr(x, "xlevels")    
      obj$call <- br_call
      obj$terms <- list(mean = xt, precision = zt, full = ft)
      obj$model <- mf
      rval$object <- obj
    }

    ## add estimating functions and model object (if desired)
    if(estfun) {
      obj$y <- y
      obj$x <- list(mean = xx, precision = xz)
      rval$estfun <- estfun.betareg(obj)
    }

    return(rval)
  }


  ## call mob
  m <- match.call(expand.dots = FALSE)
  m$formula <- mob_formula
  m$fit <- betafit
  m$control <- control
  m$link <- link
  m$link.phi <- link.phi
  m$partition <- NULL
  if("..." %in% names(m)) m[["..."]] <- NULL
  m[[1L]] <- as.call(expression(partykit::mob))[[1L]]
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  class(rval) <- c("betatree", class(rval))
  return(rval)
}


## methods
print.betatree <- function(x,
  title = "Beta regression tree", objfun = "negative log-likelihood", ...)
{
  partykit::print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.betatree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  partykit::predict.modelparty(object, newdata = newdata, type = type, ...)
}

sctest.betatree <- function(x, ...) partykit::sctest.modelparty(x, ...)

plot.betatree <- function(x, terminal_panel = partykit::node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot(partykit::as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    partykit::plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}
