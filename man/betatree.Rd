\name{betatree}
\encoding{latin1}

\alias{betatree}
\alias{coef.betatree}
\alias{logLik.betatree}
\alias{plot.betatree}
\alias{print.betatree}
\alias{sctest.betatree}
\alias{summary.betatree}
\alias{weights.betatree}

\alias{betaReg}
\alias{print.betaReg}
\alias{reweight.betaReg}
\alias{summary.betaReg}
\alias{estfun.betaReg}

\title{Beta Regression Trees}

\description{
  Fit beta regression trees via model-based recursive partitioning.
}

\usage{
betatree(formula, partition,
  data, subset = NULL, na.action = na.omit, 
  link = "logit", link.phi = "log", control = betareg.control(),
  \dots)

betaReg(control = betareg.control())
}

\arguments{
  \item{formula}{symbolic description of the model of type \code{y ~ x}
    or \code{y ~ x | z}, specifying the variables influencing mean
    and precision of \code{y}, respectively. For details see \code{\link{betareg}}.}
  \item{partition}{symbolic description of the partitioning variables,
    e.g., \code{~ p1 + p2}. The argument \code{partition} can be omitted
    if \code{formula} is a three-part formula of type \code{y ~ x | z | p1 + p2}.}  
  \item{data, subset, na.action}{arguments controlling formula processing.}
  \item{link}{character specification of the link function in
    the mean model (mu). Currently, \code{"logit"}, \code{"probit"},
    \code{"cloglog"}, \code{"cauchit"}, \code{"log"}, \code{"loglog"} are supported.
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{link.phi}{character specification of the link function in
    the precision model (phi). Currently, \code{"identity"},
    \code{"log"}, \code{"sqrt"} are supported.
    Alternatively, an object of class \code{"link-glm"} can be supplied.}
  \item{control}{a list of control arguments for the beta regression specified via
    \code{\link{betareg.control}}.}
  \item{\dots}{further control arguments for the recursive partitioning
    passed to \code{\link[party]{mob_control}}.}
}

\details{
  Beta regression trees are an application of model-based recursive partitioning
  (implemented in \code{\link[party]{mob}}, see Zeileis et al. 2008) to
  beta regression (implemented in \code{\link{betareg}}, see Cribari-Neto
  and Zeileis 2010). For plugging in \code{betareg} into \code{mob} some
  glue is required which is provided by \code{\link{betaReg}}. See also
  Grün at al. (2011) for more details.

  Various methods are provided for \code{"betatree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}, \code{summary},
  \code{coef}, etc.). The \code{plot} method employs the \code{\link[party]{node_bivplot}}
  panel-generating function.

  \code{betaReg} is a wrapper function that creates a \code{"StatModel"} object
  with certain fitting parameters passed on to \code{betareg.fit} for fitting
  beta regression models. It is the S4 interface required for plug-in to
  \code{\link[party]{mob}}. The user does not have to call this directly but
  can simply use the \code{\link{betatree}} interface.
}

\value{
  \code{betatree()} returns an object of S3 class \code{"betatree"} which is a
  list containing only a single element of S4 class \code{"mob"} (because this is
  currently not exported from the party package).

  \code{betaReg} returns an S4 object of class \code{"StatModel"} that fits beta
  regression models with the specified control arguments. When applied to data,
  it returns S3 objects of class \code{"betaReg"} inheriting from \code{"betareg"}.
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Grün, B., Kosmidis, I., and Zeileis, A. (2011).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  Working Paper 2011-22. Working Papers in Economics and Statistics,
  Research Platform Empirical and Experimental Economics, Universität Innsbruck.
  \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2011-22}

Zeileis, A., Hothorn, T., and Hornik K. (2008).
  Model-Based Recursive Partitioning.
  \emph{Journal of Computational and Graphical Statistics}, 
  \bold{17}(2), 492--514.
}

\seealso{\code{\link{betareg}}, \code{\link{betareg.fit}}, \code{\link[party]{mob}}}

\examples{
## data with two groups of dyslexic and non-dyslexic children
data("ReadingSkills", package = "betareg")
## additional random noise (not associated with reading scores)
set.seed(1071)
ReadingSkills$x1 <- rnorm(nrow(ReadingSkills))
ReadingSkills$x2 <- runif(nrow(ReadingSkills))
ReadingSkills$x3 <- factor(rnorm(nrow(ReadingSkills)) > 0)

## fit beta regression tree: in each node
##   - accurcay's mean and precision depends on iq
##   - partitioning is done by dyslexia and the noise variables x1, x2, x3
## only dyslexi is correctly selected for splitting
bt <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
  data = ReadingSkills, minsplit = 10)
plot(bt)

## inspect result
coef(bt)
sctest(bt)
summary(bt, node = 2)
summary(bt, node = 3)

## add a numerical variable with relevant information for splitting
ReadingSkills$x4 <- rnorm(nrow(ReadingSkills), c(-1.5, 1.5)[ReadingSkills$dyslexia])

bt2 <- betatree(accuracy ~ iq | iq, ~ x1 + x2 + x3 + x4,
  data = ReadingSkills, minsplit = 10)
plot(bt2)

## inspect result
coef(bt2)
sctest(bt2)
summary(bt2, node = 2)
summary(bt2, node = 3)

}

\keyword{tree}
