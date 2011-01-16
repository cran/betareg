## packages

library("betareg")
library("lmtest")
library("strucchange")


## GasolineYield

data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy_logit)

gy_loglog <- betareg(yield ~ batch + temp, data = GasolineYield,
  link = "loglog")
summary(gy_loglog)

summary(gy_logit)$pseudo.r.squared
summary(gy_loglog)$pseudo.r.squared

gy_logit2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)
summary(gy_logit2)

AIC(gy_logit, gy_logit2, gy_loglog)
lrtest(gy_logit, . ~ . + I(predict(gy_logit, type = "link")^2))
lrtest(gy_loglog, . ~ . + I(predict(gy_loglog, type = "link")^2))

gy_logit4 <- update(gy_logit, subset = -4)
coef(gy_logit, model = "precision")
coef(gy_logit4, model = "precision")

gy_loglog2 <- update(gy_loglog, link.phi = "log")
summary(gy_loglog2)$iterations


## FoodExpenditure

data("FoodExpenditure", package = "betareg")
gy_logit2 <- betareg(yield ~ batch + temp | temp, data = GasolineYield)
summary(gy_logit2)
lrtest(gy_logit, gy_logit2)

fe_lm <- lm(I(food/income) ~ income + persons, data = FoodExpenditure)
summary(fe_lm)
bptest(fe_lm)

fe_beta <- betareg(I(food/income) ~ income + persons,
  data = FoodExpenditure)
summary(fe_beta)

fe_beta2 <- betareg(I(food/income) ~ income + persons | persons,
  data = FoodExpenditure)
summary(fe_beta2)

lrtest(fe_beta, fe_beta2)
AIC(fe_beta, fe_beta2, k = log(nrow(FoodExpenditure)))

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"),
  function(x) logLik(update(fe_beta2, link = x)))


## ReadingSkills

data("ReadingSkills", package = "betareg")
rs_ols <- lm(qlogis(accuracy) ~ dyslexia * iq, data = ReadingSkills)
summary(rs_ols)

rs_beta <- betareg(accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills, hessian = TRUE)
summary(rs_beta)


## strucchange

set.seed(123)
y1 <- c(rbeta(150, 0.3 * 4, 0.7 * 4), rbeta(50, 0.5 * 4, 0.5 * 4))
y2 <- c(rbeta(100, 0.3 * 4, 0.7 * 4), rbeta(100, 0.3 * 8, 0.7 * 8))
y1_gefp <- gefp(y1 ~ 1, fit = betareg)
y2_gefp <- gefp(y2 ~ 1, fit = betareg)
sctest(y1_gefp)
sctest(y2_gefp)
