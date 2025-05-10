## -----------------------------------------------------------------------------
#| label: preliminaries
#| include: false
library("betareg")
data("GasolineYield", package = "betareg")
data("FoodExpenditure", package = "betareg")
gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")
fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)

knitr::opts_chunk$set(
  engine = "R",
  collapse = TRUE,
  comment = "##",
  message = FALSE,
  warning = FALSE,
  echo = TRUE
)
options(width = 70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE, digits = 5)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 8.5
#| fig-height: 4.5
#| out-width: 100%
#| label: fig-beta-distributions
#| fig-cap: "Probability density functions for beta distributions with varying parameters $\\mu = 0.10, 0.25, 0.50, 0.75, 0.90$ and $\\phi = 5$ (left) and $\\phi = 100$ (right)."
par(mfrow = c(1, 2), mar = c(4.1, 4.1, 4.1, 0.1))
dbeta2 <- function(x, mu, phi = 1) dbeta(x, mu * phi, (1 - mu) * phi)
x <- seq(from = 0.01, to = 0.99, length = 200)
xx <- cbind(x, x, x, x, x)

yy <- cbind(
  dbeta2(x, 0.10, 5),
  dbeta2(x, 0.25, 5),
  dbeta2(x, 0.50, 5),
  dbeta2(x, 0.75, 5),
  dbeta2(x, 0.90, 5)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "Density", main = expression(phi == 5),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.05, 12  , "0.10")
text(0.95, 12  , "0.90")
text(0.22,  2.8, "0.25")
text(0.78,  2.8, "0.75")
text(0.50,  2.3, "0.50")

yy <- cbind(
  dbeta2(x, 0.10, 100),
  dbeta2(x, 0.25, 100),
  dbeta2(x, 0.50, 100),
  dbeta2(x, 0.75, 100),
  dbeta2(x, 0.90, 100)
)
matplot(xx, yy, type = "l", xlab = "y", ylab = "", main = expression(phi == 100),
  lty = 1, col = "black", ylim = c(0, 15))
text(0.10, 14.5, "0.10")
text(0.90, 14.5, "0.90")
text(0.25,  9.8, "0.25")
text(0.75,  9.8, "0.75")
text(0.50,  8.6, "0.50")


## ----eval=FALSE---------------------------------------------------------------
# betareg(formula, data, subset, na.action, weights, offset,
#   link = "logit", link.phi = NULL, control = betareg.control(...),
#   model = TRUE, y = TRUE, x = FALSE, ...)


## -----------------------------------------------------------------------------
#| label: GasolineYield-betareg
data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy_logit)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig-height: 5.5
#| out-width: 100%
#| label: fig-GasolineYield
#| fig-cap: "Gasoline yield data from @betareg:Prater:1956: Proportion of crude oil converted to gasoline explained by temperature (in degrees Fahrenheit) at which all gasoline has vaporized and given batch (indicated by gray level). Fitted curves correspond to beta regressions `gy_loglog` with log-log link (solid, red) and `gy_logit` with logit link (dashed, blue). Both curves were evaluated at varying temperature with the intercept for batch 6 (i.e., roughly the average intercept)."
redblue <- hcl(c(0, 260), 90, 40)
plot(yield ~ temp, data = GasolineYield, type = "n",
  ylab = "Proportion of crude oil converted to gasoline",
  xlab = "Temperature at which all gasoline has vaporized",
  main = "Prater's gasoline yield data")
points(yield ~ temp, data = GasolineYield, cex = 1.75, 
  pch = 19, col = rev(gray.colors(10))[as.numeric(batch)])
points(yield ~ temp, data = GasolineYield, cex = 1.75)
legend("topleft", as.character(1:10), title = "Batch",
  col = rev(gray.colors(10)), pch = 19, bty = "n")
legend("topleft", as.character(1:10), title = "Batch", pch = 1, bty = "n")
lines(150:500, predict(gy_logit, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[2], lwd = 2, lty = 2)
lines(150:500, predict(gy_loglog, 
  newdata = data.frame(temp = 150:500, batch = "6")),
  col = redblue[1], lwd = 2)
legend("bottomright", c("log-log", "logit"),
  col = redblue, lty = 1:2, lwd = 2, bty = "n")


## -----------------------------------------------------------------------------
#| fig-width: 8.5
#| fig-height: 10
#| out-width: 100%
#| label: fig-GasolineYield-plot
#| fig-cap: "Diagnostic plots for beta regression model `gy_logit`."
par(mfrow = c(3, 2))
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
plot(gy_logit, which = 1:4, type = "pearson")
plot(gy_logit, which = 5, type = "deviance", sub.caption = "")
plot(gy_logit, which = 1, type = "deviance", sub.caption = "")


## -----------------------------------------------------------------------------
#| label: GasolineYield-update
gy_logit4 <- update(gy_logit, subset = -4)
coef(gy_logit, model = "precision")
coef(gy_logit4, model = "precision")


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-lm
data("FoodExpenditure", package = "betareg")
fe_lm <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-bptest
library("lmtest")
bptest(fe_lm)


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-betareg
fe_beta <- betareg(I(food/income) ~ income + persons,
  data = FoodExpenditure)
summary(fe_beta)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig-height: 5.5
#| out-width: 100%
#| label: fig-FoodExpenditure
#| fig-cap: "Household food expenditure data from @betareg:Griffiths+Hill+Judge:1993: Proportion of household income spent on food explained by household income and number of persons in household (indicated by gray level). Fitted curves correspond to beta regressions `fe_beta` with fixed dispersion (long-dashed, blue), `fe_beta2` with variable dispersion (solid, red), and the linear regression `fe_lin` (dashed, black). All curves were evaluated at varying income with the intercept for mean number of persons ($ = `r round(mean(FoodExpenditure$persons), digits = 2)`$)."
redblueblack <- hcl(c(0, 260, 0), c(90, 90, 0), c(40, 40, 0))
plot(I(food/income) ~ income, data = FoodExpenditure,
  xlab = "Household income", ylab = "Proportion of food expenditures",
  main = "Food expenditures data", type = "n", ylim = c(0.04, 0.57))
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75, pch = 19,
  col = rev(gray.colors(7))[persons])
points(I(food/income) ~ income, data = FoodExpenditure, cex = 1.75)
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", col = gray.colors(7), pch = 19, bty = "n")
legend("bottomleft", rev(as.character(sort(unique(FoodExpenditure$persons)))),
  title = "Persons", pch = 1, bty = "n")
lines(10:100, predict(fe_lm, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[3], lwd = 2, lty = 2)
lines(10:100, predict(fe_beta, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[2], lwd = 2, lty = 5)
lines(10:100, predict(fe_beta2, 
  newdata = data.frame(income = 10:100, persons = mean(FoodExpenditure$persons))),
  col = redblueblack[1], lwd = 2)
legend("topright", c("logit, var. disp.", "logit, fix. disp.", "lm"),
  col = redblueblack, lty = c(1, 5, 2), lwd = 2, bty = "n")


## -----------------------------------------------------------------------------
#| label: GasolineYield-phireg
gy_logit2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)


## -----------------------------------------------------------------------------
#| label: GasolineYield-phireg-coef
#| echo: false
printCoefmat(summary(gy_logit2)$coefficients$precision)


## -----------------------------------------------------------------------------
#| label: GasolineYield-lrtest
lrtest(gy_logit, gy_logit2)


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-betareg2
fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-comparison
lrtest(fe_beta, fe_beta2)
AIC(fe_beta, fe_beta2, k = log(nrow(FoodExpenditure)))


## -----------------------------------------------------------------------------
#| label: GasolineYield-loglog
gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")


## -----------------------------------------------------------------------------
#| label: GasolineYield-Rsquared
summary(gy_logit)$pseudo.r.squared
summary(gy_loglog)$pseudo.r.squared


## -----------------------------------------------------------------------------
#| label: GasolineYield-AIC
AIC(gy_logit, gy_logit2, gy_loglog)


## -----------------------------------------------------------------------------
#| label: GasolineYield-reset
lrtest(gy_logit, . ~ . + I(predict(gy_logit, type = "link")^2))
lrtest(gy_loglog, . ~ . + I(predict(gy_loglog, type = "link")^2))


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig-height: 5.5
#| out-width: 100%
#| label: fig-GasolineYield-diagnostics
#| fig-cap: "Scatterplot comparing the absolute raw residuals from beta regression modes with log-log link (x-axis) and logit link (y-axis)."
plot(abs(residuals(gy_loglog, type = "response")),
  abs(residuals(gy_logit, type = "response")))
abline(0, 1, lty = 2)


## -----------------------------------------------------------------------------
#| label: GasolineYield-loglog2
gy_loglog2 <- update(gy_loglog, link.phi = "log")
summary(gy_loglog2)$iterations


## -----------------------------------------------------------------------------
#| label: FoodExpenditure-links
sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
  function(x) logLik(update(fe_beta2, link = x)))


## -----------------------------------------------------------------------------
#| label: ReadingSkills-eda
#| echo: false
#| results: hide
data("ReadingSkills", package = "betareg")
rs_accuracy <- format(round(with(ReadingSkills, tapply(accuracy, dyslexia, mean)), digits = 3))


## -----------------------------------------------------------------------------
#| label: ReadingSkills-ols
data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
coeftest(rs_ols)


## -----------------------------------------------------------------------------
#| label: ReadingSkills-beta
rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
coeftest(rs_beta)


## -----------------------------------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig.height: 5.5
#| out-width: 100%
#| label: fig-ReadingSkills
#| fig-cap: "Reading skills data from @betareg:Smithson+Verkuilen:2006 : Linearly transformed reading accuracy by IQ score and dyslexia status (control, blue vs. dyslexic, red). Fitted curves correspond to beta regression `rs_beta` (solid) and OLS regression with logit-transformed dependent variable `rs_ols` (dashed)."
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
#| label: strucchange-data
suppressWarnings(RNGversion("3.5.0"))
set.seed(123)
y1 <- c(rbeta(150, 0.3 * 4, 0.7 * 4), rbeta(50, 0.5 * 4, 0.5 * 4))
y2 <- c(rbeta(100, 0.3 * 4, 0.7 * 4), rbeta(100, 0.3 * 8, 0.7 * 8))


## -----------------------------------------------------------------------------
#| label: strucchange-gefp
library("strucchange")
y1_gefp <- gefp(y1 ~ 1, fit = betareg)
y2_gefp <- gefp(y2 ~ 1, fit = betareg)


## -----------------------------------------------------------------------------
#| fig-width: 6.5
#| fig-height: 6
#| label: fig-strucchange1
#| fig-cap: "Structural change tests for artificial data `y1` with change in $\\mu$."
plot(y1_gefp, aggregate = FALSE)


## -----------------------------------------------------------------------------
#| fig-width: 6.5
#| fig-height: 6
#| label: fig-strucchange2
#| fig-cap: "Structural change tests for artificial data `y2` with change in $\\phi$."
plot(y2_gefp, aggregate = FALSE)

