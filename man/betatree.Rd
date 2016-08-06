\name{betatree}
\encoding{latin1}

\alias{betatree}
\alias{plot.betatree}
\alias{print.betatree}
\alias{predict.betatree}
\alias{sctest.betatree}

\title{Beta Regression Trees}

\description{
  Fit beta regression trees via model-based recursive partitioning.
}

\usage{
betatree(formula, partition,
  data, subset = NULL, na.action = na.omit, weights, offset, cluster,
  link = "logit", link.phi = "log", control = betareg.control(),
  \dots)
}

\arguments{
  \item{formula}{symbolic description of the model of type \code{y ~ x}
    or \code{y ~ x | z}, specifying the variables influencing mean
    and precision of \code{y}, respectively. For details see \code{\link{betareg}}.}
  \item{partition}{symbolic description of the partitioning variables,
    e.g., \code{~ p1 + p2}. The argument \code{partition} can be omitted
    if \code{formula} is a three-part formula of type \code{y ~ x | z | p1 + p2}.}  
  \item{data, subset, na.action, weights, offset, cluster}{arguments controlling
    data/model processing passed to \code{\link[partykit]{mob}}.}
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
    passed to \code{\link[partykit]{mob_control}}.}
}

\details{
  Beta regression trees are an application of model-based recursive partitioning
  (implemented in \code{\link[partykit]{mob}}, see Zeileis et al. 2008) to
  beta regression (implemented in \code{\link{betareg}}, see Cribari-Neto
  and Zeileis 2010). See also Grün at al. (2012) for more details.

  Various methods are provided for \code{"betatree"} objects, most of them
  inherit their behavior from \code{"mob"} objects (e.g., \code{print}, \code{summary},
  \code{coef}, etc.). The \code{plot} method employs the \code{\link[partykit]{node_bivplot}}
  panel-generating function.
}

\value{
  \code{betatree()} returns an object of S3 class \code{"betatree"} which
  inherits from \code{"modelparty"}.
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Grün, B., Kosmidis, I., and Zeileis, A. (2012).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  \emph{Journal of Statistical Software}, \bold{48}(11), 1--25.
  \url{http://www.jstatsoft.org/v48/i11/}.

Zeileis, A., Hothorn, T., and Hornik K. (2008).
  Model-Based Recursive Partitioning.
  \emph{Journal of Computational and Graphical Statistics}, 
  \bold{17}(2), 492--514.
}

\seealso{\code{\link{betareg}}, \code{\link{betareg.fit}}, \code{\link[partykit]{mob}}}

\examples{
options(digits = 4)

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
## only dyslexia is correctly selected for splitting
bt <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
  data = ReadingSkills, minsize = 10)
plot(bt)

## inspect result
coef(bt)
summary(bt, node = 2)
summary(bt, node = 3)
library("strucchange")
sctest(bt)

## add a numerical variable with relevant information for splitting
ReadingSkills$x4 <- rnorm(nrow(ReadingSkills), c(-1.5, 1.5)[ReadingSkills$dyslexia])

bt2 <- betatree(accuracy ~ iq | iq, ~ x1 + x2 + x3 + x4,
  data = ReadingSkills, minsize = 10)
plot(bt2)

## inspect result
coef(bt2)
sctest(bt2)
summary(bt2, node = 2)
summary(bt2, node = 3)
}

\keyword{tree}
