\name{ReadingSkills}
\encoding{latin1}

\alias{ReadingSkills}

\title{Dyslexia and IQ Predicting Reading Accuracy}

\description{
Data for assessing the contribution of non-verbal IQ to
children's reading skills in dyslexic and non-dyslexic children.
}

\usage{data("ReadingSkills")}

\format{
A data frame containing 44 observations on 3 variables.
\describe{
  \item{accuracy}{reading score scaled to the open unit interval (see below).}
  \item{dyslexia}{factor. Is the child dyslexic? (A sum contrast rather
    than treatment contrast is employed.)}
  \item{iq}{non-verbal intelligence quotient transformed to z-scores.}
}
}

\details{
  The data were collected by Pammer and Kevan (2004) and employed by
  Smithson and Verkuilen (2006). The original reading accuracy score was transformed
  by Smithson and Verkuilen (2006) so that \code{accuracy} is in the open unit
  interval (0, 1) and beta regression can be employed. First, the original accuracy
  was scaled using the minimal and maximal score (\code{a} and \code{b}, respectively)
  that can be obtained in the test: \code{(original_accuracy - a) / (b - a)}
  (\code{a} and \code{b} are not provided). Subsequently, the scaled score is transformed
  to the unit interval using a continuity correction: \code{(scaled_accuracy * (n-1) - 0.5) / n}
  (either with some rounding or using \code{n = 50} rather than 44).
}

\source{
Example 3 from
\url{http://psychology3.anu.edu.au/people/smithson/details/betareg/betareg.html}
}

\references{
Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \url{http://www.jstatsoft.org/v34/i02/}.

Grün, B., Kosmidis, I., and Zeileis, A. (2011).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  Working Paper 2011-22. Working Papers in Economics and Statistics,
  Research Platform Empirical and Experimental Economics, Universitaet Innsbruck.
  \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2011-22}

Pammer, K., and Kevan, A. (2004).
  The Contribution of Visual Sensitivity, Phonological Processing
  and Non-Verbal IQ to Children's Reading.
  \emph{Unpublished manuscript}, The Australian National University, Canberra.

Smithson, M., and Verkuilen, J. (2006).
  A Better Lemon Squeezer? Maximum-Likelihood Regression with
  Beta-Distributed Dependent Variables.
  \emph{Psychological Methods}, \bold{11}(7), 54--71.
}

\seealso{\code{\link{betareg}}, \code{\link{MockJurors}}, \code{\link{StressAnxiety}}}

\examples{
data("ReadingSkills", package = "betareg")

## Smithson & Verkuilen (2006, Table 5)
## OLS regression
## (Note: typo in iq coefficient: 0.3954 instead of 0.3594)
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
summary(rs_ols)
## Beta regression (with numerical rather than analytic standard errors)
## (Note: Smithson & Verkuilen erroneously compute one-sided p-values)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
summary(rs_beta)

## visualization
plot(accuracy ~ iq, data = ReadingSkills, col = as.numeric(dyslexia), pch = 19)
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd))
lines(nd$iq, plogis(predict(rs_ols, nd)), lty = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = 2, lty = 2)

## see demo("SmithsonVerkuilen2006", package = "betareg") for more details
}

\keyword{datasets}
