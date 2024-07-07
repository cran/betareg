\name{StressAnxiety}
\alias{StressAnxiety}

\title{Dependency of Anxiety on Stress}

\description{
Stress and anxiety among nonclinical women in Townsville, Queensland, Australia.
}

\usage{data("StressAnxiety", package = "betareg")}

\format{
A data frame containing 166 observations on 2 variables.
\describe{
  \item{stress}{score, linearly transformed to the open unit
    interval (see below).}
  \item{anxiety}{score, linearly transformed to the open unit
    interval (see below).}
}
}

\details{
  Both variables were assess on the Depression Anxiety Stress Scales, ranging
  from 0 to 42. Smithson and Verkuilen (2006) transformed these to the open
  unit interval (without providing details about this transformation).
}

\source{
Example 2 from Smithson and Verkuilen (2006) supplements.
}

\references{
Smithson, M., and Verkuilen, J. (2006).
  A Better Lemon Squeezer? Maximum-Likelihood Regression with
  Beta-Distributed Dependent Variables.
  \emph{Psychological Methods}, \bold{11}(7), 54--71.
}

\seealso{\code{\link{betareg}}, \code{\link{MockJurors}}, \code{\link{ReadingSkills}}}

\examples{
data("StressAnxiety", package = "betareg")
StressAnxiety <- StressAnxiety[order(StressAnxiety$stress),]

## Smithson & Verkuilen (2006, Table 4)
sa_null <- betareg(anxiety ~ 1 | 1,
  data = StressAnxiety, hessian = TRUE)
sa_stress <- betareg(anxiety ~ stress | stress,
  data = StressAnxiety, hessian = TRUE)
summary(sa_null)
summary(sa_stress)
AIC(sa_null, sa_stress)
1 - as.vector(logLik(sa_null)/logLik(sa_stress))

## visualization
attach(StressAnxiety)
plot(jitter(anxiety) ~ jitter(stress),
  xlab = "Stress", ylab = "Anxiety",
  xlim = c(0, 1), ylim = c(0, 1))
lines(lowess(anxiety ~ stress))
lines(fitted(sa_stress) ~ stress, lty = 2)
lines(fitted(lm(anxiety ~ stress)) ~ stress, lty = 3)
legend("topleft", c("lowess", "betareg", "lm"), lty = 1:3, bty = "n")
detach(StressAnxiety)

## see demo("SmithsonVerkuilen2006", package = "betareg") for more details
}

\keyword{datasets}
