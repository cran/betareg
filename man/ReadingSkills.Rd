\name{ReadingSkills}
\encoding{utf8}

\alias{ReadingSkills}

\title{Dyslexia and IQ Predicting Reading Accuracy}

\description{
Data for assessing the contribution of non-verbal IQ to
children's reading skills in dyslexic and non-dyslexic children.
}

\usage{data("ReadingSkills", package = "betareg")}

\format{
A data frame containing 44 observations on 3 variables.
\describe{
  \item{accuracy}{numeric. Reading score with maximum restricted to be
    0.99 rather than 1 (see below).}
  \item{dyslexia}{factor. Is the child dyslexic? (A sum contrast rather
    than treatment contrast is employed.)}
  \item{iq}{numeric. Non-verbal intelligence quotient transformed to z-scores.}
  \item{accuracy1}{numeric. Unrestricted reading score with a maximum of 1
    (see below).}
}
}

\details{
  The data were collected by Pammer and Kevan (2004) and employed by
  Smithson and Verkuilen (2006). The original reading accuracy score was transformed
  by Smithson and Verkuilen (2006) so that \code{accuracy} is in the open unit
  interval (0, 1) and beta regression can be employed. First, the original accuracy
  was scaled using the minimal and maximal score (\code{a} and \code{b}, respectively)
  that can be obtained in the test: \code{accuracy1 = (original_accuracy - a) / (b - a)}
  (\code{a} and \code{b} are not provided). Subsequently, \code{accuracy} was obtained
  from \code{accuracy1} by replacing all observations with a value of 1 with 0.99.
  
  Kosmidis and Zeileis (2024) propose to investigate the original unrestricted
  \code{accuracy1} variable using their extended-support beta mixture regression.
}

\source{
Example 3 from Smithson and Verkuilen (2006) supplements.
}

\references{
Cribari-Neto F, Zeileis A (2010). Beta Regression in R.
  \emph{Journal of Statistical Software}, \bold{34}(2), 1--24.
  \doi{10.18637/jss.v034.i02}

Grün B, Kosmidis I, Zeileis A (2012).
  Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
  \emph{Journal of Statistical Software}, \bold{48}(11), 1--25.
  \doi{10.18637/jss.v048.i11}

Kosmidis I, Zeileis A (2024).
  Extended-Support Beta Regression for [0, 1] Responses.
  2409.07233, \emph{arXiv.org E-Print Archive}.
  \doi{10.48550/arXiv.2409.07233}

Pammer K, Kevan A (2004).
  The Contribution of Visual Sensitivity, Phonological Processing
  and Non-Verbal IQ to Children's Reading.
  \emph{Unpublished manuscript}, The Australian National University, Canberra.

Smithson M, Verkuilen J (2006).
  A Better Lemon Squeezer? Maximum-Likelihood Regression with
  Beta-Distributed Dependent Variables.
  \emph{Psychological Methods}, \bold{11}(7), 54--71.
}

\seealso{\code{\link{betareg}}, \code{\link{MockJurors}}, \code{\link{StressAnxiety}}}

\examples{
options(digits = 4)
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

## Extended-support beta mixture regression (Kosmidis & Zeileis 2024)
rs_xbx <- betareg(accuracy1 ~ dyslexia * iq | dyslexia + iq, data = ReadingSkills)
summary(rs_xbx)

## Coefficients in XBX are typically somewhat shrunken compared to beta
cbind(XBX = coef(rs_xbx), Beta = c(coef(rs_beta), NA))

## Visualization
plot(accuracy1 ~ iq, data = ReadingSkills, col = c(4, 2)[dyslexia], pch = 19)
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_xbx, nd), col = 4)
lines(nd$iq, predict(rs_beta, nd), col = 4, lty = 5)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = 4, lty = 3)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_xbx, nd), col = 2)
lines(nd$iq, predict(rs_beta, nd), col = 2, lty = 5)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = 2, lty = 3)
legend("topleft", c("Dyslexia: no", "Dyslexia: yes", "OLS", "XBX", "Beta"),
  lty = c(0, 0, 3, 1, 5), pch = c(19, 19, NA, NA, NA), col = c(4, 2, 1, 1, 1), bty = "n")

## see demo("SmithsonVerkuilen2006", package = "betareg") for further details
}

\keyword{datasets}
