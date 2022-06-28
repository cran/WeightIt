## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warnings = FALSE, messages = FALSE)
set.seed(1000)

## -----------------------------------------------------------------------------
data("lalonde", package = "cobalt")
head(lalonde)

## -----------------------------------------------------------------------------
library("cobalt")
bal.tab(treat ~ age + educ + race + married + nodegree + re74 + re75,
        data = lalonde, estimand = "ATT", thresholds = c(m = .05))

## -----------------------------------------------------------------------------
library("WeightIt")
W.out <- weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                  data = lalonde, estimand = "ATT", method = "ps")
W.out #print the output

## -----------------------------------------------------------------------------
summary(W.out)

## -----------------------------------------------------------------------------
bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

## -----------------------------------------------------------------------------
W.out <- weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                  data = lalonde, estimand = "ATT", method = "ebal")
summary(W.out)

## -----------------------------------------------------------------------------
bal.tab(W.out, stats = c("m", "v"), thresholds = c(m = .05))

## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(eval = requireNamespace("survey", quietly = TRUE))

## ---- message=FALSE-----------------------------------------------------------
library(survey)
d.w <- svydesign(~1, weights = W.out$weights, data = lalonde)
fit <- svyglm(re78 ~ treat, design = d.w)
coef(fit)

## -----------------------------------------------------------------------------
#Robust standard errors and confidence intervals
summary(fit)
confint(fit)

## ---- include=FALSE, eval=TRUE------------------------------------------------
knitr::opts_chunk$set(eval = requireNamespace("boot", quietly = TRUE))

## ---- warning=FALSE, message=FALSE--------------------------------------------
#Bootstrapping
library("boot")
est.fun <- function(data, index) {
  W.out <- weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                    data = data[index,], estimand = "ATT", method = "ebal")
  fit <- glm(re78 ~ treat, data = data[index,], weights = W.out$weights)
  return(coef(fit)["treat"])
}
boot.out <- boot(est.fun, data = lalonde, R = 999)
boot.ci(boot.out, type = "bca") #type shouldn't matter so much

## ---- include=FALSE, eval=TRUE------------------------------------------------
knitr::opts_chunk$set(eval = all(sapply(c("twang", "survey"), requireNamespace, quietly = TRUE)))

## -----------------------------------------------------------------------------
data("iptwExWide", package = "twang")
head(iptwExWide)

## -----------------------------------------------------------------------------
library("cobalt") #if not already attached
bal.tab(list(tx1 ~ age + gender + use0,
             tx2 ~ tx1 + use1 + age + gender + use0,
             tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
        data = iptwExWide, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .all)

## -----------------------------------------------------------------------------
Wmsm.out <- weightitMSM(list(tx1 ~ age + gender + use0,
                             tx2 ~ tx1 + use1 + age + gender + use0,
                             tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
                        data = iptwExWide, method = "ps",
                        stabilize = TRUE)
Wmsm.out

## -----------------------------------------------------------------------------
summary(Wmsm.out)

## -----------------------------------------------------------------------------
bal.tab(Wmsm.out, stats = c("m", "ks"), thresholds = c(m = .05),
        which.time = .none)

## ---- message=FALSE-----------------------------------------------------------
library("survey")
d.w.msm <- svydesign(~1, weights = Wmsm.out$weights,
                     data = iptwExWide)
full.fit <- svyglm(outcome ~ tx1*tx2*tx3, design = d.w.msm)
main.effects.fit <- svyglm(outcome ~ tx1 + tx2 + tx3, design = d.w.msm)
anova(full.fit, main.effects.fit)

## -----------------------------------------------------------------------------
cum.fit <- svyglm(outcome ~ I(tx1+tx2+tx3), design = d.w.msm)
anova(main.effects.fit, cum.fit)
anova(full.fit, cum.fit)

## -----------------------------------------------------------------------------
summary(cum.fit)
confint(cum.fit)

## ---- include=FALSE, eval=TRUE------------------------------------------------
knitr::opts_chunk$set(eval = TRUE)

