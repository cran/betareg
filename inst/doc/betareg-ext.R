## -----------------------------------------------------------------------------
#| label: preliminaries
#| include: false
library("betareg")

knitr::opts_chunk$set(
  engine = "R",
  collapse = TRUE,
  comment = "##",
  message = FALSE,
  warning = FALSE,
  echo = TRUE
)
options(width = 70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE, digits = 5)

combine <- function(x, sep, width) {
  cs <- cumsum(nchar(x))
  remaining <- if (any(cs[-1] > width)) combine(x[c(FALSE, cs[-1] > width)], sep, width)
  c(paste(x[c(TRUE, cs[-1] <= width)], collapse= sep), remaining)
}
prettyPrint <- function(x, sep = " ", linebreak = "\n\t", width = getOption("width")) {
  x <- strsplit(x, sep)[[1]]
  paste(combine(x, sep, width), collapse = paste(sep, linebreak, collapse = ""))
}
cache <- FALSE
enumerate <- function(x) paste(paste(x[-length(x)], collapse = ", "), x[length(x)], sep = " and ")
betamix_methods <- enumerate(paste("`", gsub("\\.betamix", "", as.character(methods(class = "betamix"))), "`", sep = ""))


## ----eval=FALSE---------------------------------------------------------------
# betareg(formula, data, subset, na.action, weights, offset,
#   link = "logit", link.phi = NULL, type = c("ML", "BC", "BR"),
#   control = betareg.control(...), model = TRUE, y = TRUE, x = FALSE, ...)


## ----eval=FALSE---------------------------------------------------------------
# betatree(formula, partition, data, subset, na.action, weights, offset,
#   link = "logit", link.phi = "log", control = betareg.control(), ...)


## ----eval=FALSE---------------------------------------------------------------
# betamix(formula, data, k, fixed, subset, na.action,
#   link = "logit", link.phi = "log", control = betareg.control(...),
#   FLXconcomitant = NULL, extra_components,
#   verbose = FALSE, ID, nstart = 3, FLXcontrol = list(), cluster = NULL,
#   which = "BIC", ...)


## ----include=FALSE------------------------------------------------------------
data("ReadingSkills", package = "betareg")
mean_accuracy <-
  format(round(with(ReadingSkills, tapply(accuracy, dyslexia, mean)), digits = 3),
         nsmall = 3)
mean_iq <-
  format(round(with(ReadingSkills, tapply(iq, dyslexia, mean)), digits = 3),
         nsmall = 3)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig-height: 5.5
#| out-width: 100%
#| label: fig-ReadingSkills
#| fig-cap: "Reading skills data from @betareg:Smithson+Verkuilen:2006. Linearly transformed reading accuracy by IQ score and dyslexia status (control, blue vs. dyslexic, red). Fitted curves correspond to beta regression (solid) and OLS regression with logit-transformed dependent variable (dashed)."
data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
cl1 <- hcl(c(260, 0), 90, 40)
cl2 <- hcl(c(260, 0), 10, 95)
plot(accuracy ~ iq, data = ReadingSkills, col = cl2[as.numeric(dyslexia)],
  main = "Reading skills data", xlab = "IQ score", ylab = "Reading accuracy",
  pch = c(19, 17)[as.numeric(dyslexia)], cex = 1.5)
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5,
  pch = (1:2)[as.numeric(dyslexia)], col = cl1[as.numeric(dyslexia)])
nd <- data.frame(dyslexia = "no", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[1], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[1], lty = 2, lwd = 2)
nd <- data.frame(dyslexia = "yes", iq = -30:30/10)
lines(nd$iq, predict(rs_beta, nd), col = cl1[2], lwd = 2)
lines(nd$iq, plogis(predict(rs_ols, nd)), col = cl1[2], lty = 2, lwd = 2)
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(19, 17, NA, NA), lwd = 2,
  col = c(cl2, 1, 1), bty = "n")
legend("topleft", c("control", "dyslexic", "betareg", "lm"),
  lty = c(NA, NA, 1:2), pch = c(1, 2, NA, NA),
  col = c(cl1, NA, NA), bty = "n")


## -----------------------------------------------------------------------------
#| label: ReadingSkills-bias
data("ReadingSkills", package = "betareg")
rs_f <- accuracy ~ dyslexia * iq | dyslexia * iq
rs_ml <- betareg(rs_f, data = ReadingSkills, type = "ML")
rs_bc <- betareg(rs_f, data = ReadingSkills, type = "BC")
rs_br <- betareg(rs_f, data = ReadingSkills, type = "BR")


## -----------------------------------------------------------------------------
#| label: ReadingSkills-bias-table
#| echo: false
#| results: asis
rs_list <- list(rs_ml, rs_bc, rs_br)
cf <- paste("$", sapply(round(sapply(rs_list, coef), digits = 3), format, nsmall = 3), "$", sep = "")
se <- paste("$", format(round(sapply(rs_list, function(x) sqrt(diag(vcov(x)))), digits = 3), nsmall = 3), "$", sep = "")
ll <- paste("$", format(round(sapply(rs_list, logLik), digits = 3), nsmall = 3), "$", sep = "")
cfse <- matrix(as.vector(rbind(cf, se)), ncol = 3)
cfse <- cbind(
  c("Mean", rep("", 7), "Precision", rep("", 7)),
  rep(as.vector(rbind(c("(Intercept)", "`dyslexia`", "`iq`", "`dyslexia:iq`"), "")), 2),
  cfse)
cfse <- rbind(cfse, c("Log-likelihood", "", ll))
knitr::kable(cfse, align = c("l", "l", "r", "r", "r"), col.names = c("", "", "Maximum likelihood", "Bias correction", "Bias reduction"))


## -----------------------------------------------------------------------------
#| echo: false
#| fig-height: 6.5
#| fig-width: 6.5
#| out-width: 100%
#| label: fig-readingskillsbias
#| fig-cap: "Scatterplots of the logarithm of the estimated precision parameters $\\log(\\phi_i)$ based on the maximum likelihood, bias-corrected and bias-reduced estimates. The dashed black line is the main diagonal, the solid red line is a scatterplot smoother."
pr_phi <- sapply(list("Maximum likelihood" = rs_ml,
                      "Bias correction" = rs_bc,
                      "Bias reduction" = rs_br), predict, type = "precision")
pairs(log(pr_phi), panel = function(x, y, ...) {
  panel.smooth(x, y, ...)
  abline(0, 1, lty = 2)
  })


## -----------------------------------------------------------------------------
#| label: ReadingSkills-noise
#| echo: true
suppressWarnings(RNGversion("3.5.0"))
set.seed(1071)
n <- nrow(ReadingSkills)
ReadingSkills$x1 <- rnorm(n)
ReadingSkills$x2 <- runif(n)
ReadingSkills$x3 <- factor(sample(0:1, n, replace = TRUE))


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree
rs_tree <- betatree(accuracy ~ iq | iq, ~ dyslexia + x1 + x2 + x3,
  data = ReadingSkills, minsize = 10)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree2
#| echo: true
#| eval: false
# rs_tree <- betatree(accuracy ~ iq | iq | dyslexia + x1 + x2 + x3,
#   data = ReadingSkills, minsize = 10)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree3
#| eval: false
# plot(rs_tree)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-height: 7
#| fig-width: 10
#| out-width: 100%
#| label: fig-betatree
#| fig-cap: "Partitioned beta regression model for the `ReadingSkills` data."
plot(rs_tree)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree-coef
coef(rs_tree)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree4
#| echo: true
rs_tree


## -----------------------------------------------------------------------------
#| label: ReadingSkills-tree-sctest
library("strucchange")
sctest(rs_tree)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-mix
rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
  extra_components = extraComponent(type = "uniform",
    coef = 0.99, delta = 0.01), nstart = 10)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-mix3
rs_mix


## -----------------------------------------------------------------------------
#| label: ReadingSkills-mix4
summary(rs_mix)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-height: 5.5
#| fig-width: 10
#| out-width: 100%
#| label: fig-betamix
#| fig-cap: "Fitted regression lines for the mixture model with three components and the observations shaded according to their posterior probabilities (left). Fitted regression lines for the partitioned beta regression model with shading according to the observed `dyslexic` variable where nondyslexic and dyslexic children are in blue and red, respectively (right)."
par(mfrow = c(1, 2))
ix <- as.numeric(ReadingSkills$dyslexia)
prob <- 2 * (posterior(rs_mix)[cbind(seq_along(ix), clusters(rs_mix))] - 0.5)
col3 <- hcl(c(0, 260, 130), 65, 45, fixup = FALSE)
col1 <- col3[clusters(rs_mix)]
col2 <- hcl(c(0, 260, 130)[clusters(rs_mix)], 65 * abs(prob)^1.5, 95 - 50 * abs(prob)^1.5, fixup = FALSE)
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19, cex = 1.5,
  xlim = c(-2, 2), main = "Mixture model (dyslexia unobserved)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)
iq <- -30:30/10
cf <- rbind(coef(rs_mix, model = "mean", component = 1:2), c(qlogis(0.99), 0))
for(i in 1:3) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])

ix <- as.numeric(ReadingSkills$dyslexia)
col1 <- hcl(c(260, 0), 90, 40)[ix]
col2 <- hcl(c(260, 0), 10, 95)[ix]
plot(accuracy ~ iq, data = ReadingSkills, col = col2, pch = 19,
  cex = 1.5, xlim = c(-2, 2), main = "Partitioned model (dyslexia observed)")
points(accuracy ~ iq, data = ReadingSkills, cex = 1.5, pch = 1, col = col1)

cf <- coef(rs_tree, model = "mean")
col3 <- hcl(c(260, 0), 90, 40)
for(i in 1:2) lines(iq, plogis(cf[i, 1] + cf[i, 2] * iq), lwd = 2, col = col3[i])


## -----------------------------------------------------------------------------
#| label: ReadingSkills-mix5
table(clusters(rs_mix), ReadingSkills$dyslexia)


## -----------------------------------------------------------------------------
#| label: GasolineYield-bias
data("GasolineYield", package = "betareg")
gy <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp, data = GasolineYield, type = x))


## -----------------------------------------------------------------------------
#| label: GasolineYield-bias-table
#| echo: false
#| results: asis
cf <- matrix(paste("$", sapply(round(sapply(gy, coef), digits = 5), format, nsmall = 5), "$", sep = ""), ncol = 3)
se <- matrix(gsub(" ", "",
  paste("$", format(round(sapply(gy, function(x) sqrt(diag(vcov(x)))), digits = 5), nsmall = 5), "$", sep = ""),
  fixed = TRUE), ncol = 3)
cfse <- cbind(c(paste("$\\beta_{", 1:11, "}$", sep = ""), "$\\phi$"), cf[,1], se[,1], cf[,2], se[,2], cf[,3], se[,3])
knitr::kable(cfse, align = c("l", rep("r", 6)), col.names = c("", "Maximum likelihood", "", "Bias correction", "", "Bias reduction", ""))


## -----------------------------------------------------------------------------
#| label: GasolineYield-phi
sapply(gy, coef, model = "precision")


## -----------------------------------------------------------------------------
#| label: GasolineYield-phi-loglik
sapply(gy, logLik)


## -----------------------------------------------------------------------------
#| label: GasolineYield-bias2
data("GasolineYield", package = "betareg")
gy2 <- lapply(c("ML", "BC", "BR"), function(x)
  betareg(yield ~ batch + temp | 1, data = GasolineYield, type = x))
sapply(gy2, logLik)


## -----------------------------------------------------------------------------
#| label: GasolineYield-bias-table2
#| echo: false
#| results: asis
cf <- matrix(paste("$", sapply(round(sapply(gy2, coef), digits = 5), format, nsmall = 5), "$", sep = ""), ncol = 3)
se <- matrix(gsub(" ", "",
  paste("$", format(round(sapply(gy2, function(x) sqrt(diag(vcov(x)))), digits = 5), nsmall = 5), "$", sep = ""),
  fixed = TRUE), ncol = 3)
cfse <- cbind(c(paste("$\\beta_{", 1:11, "}$", sep = ""), "$\\log\\phi$"), cf[,1], se[,1], cf[,2], se[,2], cf[,3], se[,3])
knitr::kable(cfse, align = c("l", rep("r", 6)), col.names = c("", "Maximum likelihood", "", "Bias correction", "", "Bias reduction", ""))

