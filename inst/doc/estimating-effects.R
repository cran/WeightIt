## ----include = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=T)
options(width = 200, digits= 4)

me_ok <- requireNamespace("marginaleffects", quietly = TRUE)
su_ok <- requireNamespace("survival", quietly = TRUE)

## ----include = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generating data similar to Austin (2009) for demonstrating treatment effect estimation
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

gen_Ac <- function(X) {
  LP_A <- -1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
  LP_A + rlogis(nrow(X))
}

#~20% treated
gen_A <- function(Ac) {
  1 * (Ac > 0)
}

gen_Am <- function(A) {
  factor(ifelse(A == 1, "T", sample(c("C1", "C2"), length(A), TRUE)))
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + 
    .5*A*X[,1] + .5*A*X[,2] - .25*A*X[,3] + A*(X[,5] - .5) +
    rnorm(length(A), 0, 5)
}
#Conditional:
#  MD: 2
#Marginal:
#  MD: 2

# Binary outcome
gen_Y_B <- function(A, X) {
  LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  P_B <- plogis(LP_B)
  rbinom(length(A), 1, P_B)
}
#Conditional:
#  OR:   2.4
#  logOR: .875
#Marginal:
#  RD:    .144
#  RR:   1.54
#  logRR: .433
#  OR:   1.92
#  logOR  .655

# Survival outcome
gen_Y_S <- function(A, X) {
  LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
}
#Conditional:
#  HR:   2.4
#  logHR: .875
#Marginal:
#  HR:   1.57
#  logHR: .452

set.seed(19599)

n <- 2000
X <- gen_X(n)
Ac <- gen_Ac(X)
A <- gen_A(Ac)
Am <- gen_Am(A)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, Am, Ac, X, Y_C, Y_B, Y_S)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(d)

## ----message=FALSE,warning=FALSE, include=FALSE, eval=!me_ok------------------------------------------------------------------------------------------------------------------------------------------
#  library("WeightIt")

## ----message=FALSE,warning=FALSE, eval=me_ok----------------------------------------------------------------------------------------------------------------------------------------------------------
library("WeightIt")
library("marginaleffects")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PS weighting for the ATE with a logistic regression PS
W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "glm", estimand = "ATE")
W

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Linear model with covariates
fit <- lm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
                    data = d, weightit = W)

## ----eval=me_ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_comparisons(fit, variables = "A")

## ----eval=me_ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_predictions(fit, variables = "A")

## ----eval=me_ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Logistic regression model with covariates
fit <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
                    data = d, weightit = W,
                    family = binomial)

#Compute effects; RR and confidence interval
avg_comparisons(fit,
                variables = "A",
                comparison = "lnratioavg",
                transform = "exp")

## ----eval=su_ok, message=F----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cox Regression for marginal HR
fit <- coxph_weightit(survival::Surv(Y_S) ~ A, data = d,
                      weightit = W)

#Log HR estimates
summary(fit)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#HR and CIs; requested by exponentiating log HRs
summary(fit, ci = TRUE, transform = "exp")

## ----eval= FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Estimate the balancing weights, with sampling weights called "sw"
#  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 +
#                  X6 + X7 + X8 + X9, data = d,
#                method = "glm", estimand = "ATE",
#                s.weights = "sw")
#  
#  #Fit the outcome model, with clustering variable "clu"
#  fit <- glm_weightit(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 +
#                             X6 + X7 + X8 + X9),
#                      data = d, weightit = W,
#                      cluster = ~clu)
#  
#  #Compute the ATE, include sampling weights in the estimation
#  avg_comparisons(fit,
#                  variables = "A",
#                  wts = W$s.weights)

## ----eval= FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  library(survey)
#  
#  #Multiply sampling weights and estimated weights
#  d$weights <- W$weights * d$sw
#  
#  #Declare a survey design using the combined weights with
#  #appropriate clustering
#  des <- svydesign(~clu, weights = ~weights, data = d)
#  
#  #Fit the outcome model
#  fit <- svyglm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 +
#                             X6 + X7 + X8 + X9),
#                design = des)
#  
#  #G-computation for the difference in means, including sampling weights
#  avg_comparisons(fit,
#                  variables = "A",
#                  wts = d$sw)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(d$Am)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "ebal", estimand = "ATT",
              focal = "T")
W

## ----eval=me_ok---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fit the outcome model
fit <- lm_weightit(Y_C ~ Am * (X1 + X2 + X3 + X4 + X5 + 
                                 X6 + X7 + X8 + X9),
                   data = d, weightit = W)

#G-computation
p <- avg_predictions(fit,
                     variables = "Am",
                     newdata = subset(Am == "T"))
p

hypotheses(p, "revpairwise")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
W <- weightit(Ac ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              moments = 2, int = TRUE,
              method = "ebal")
W

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fit the outcome model
fit <- lm_weightit(Y_C ~ splines::ns(Ac, df = 4) *
                     (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
                   data = d, weightit = W)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Represenative values of Ac:
values <- with(d, seq(quantile(Ac, .1),
                      quantile(Ac, .9),
                      length.out = 31))

#G-computation
p <- avg_predictions(fit,
                     variables = list(Ac = values))

## ----fig.height=3.5, fig.width=7----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("ggplot2")
ggplot(p, aes(x = Ac)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  labs(x = "Ac", y = "E[Y|A]") +
  theme_bw()

## ----fig.height=3.5, fig.width=7----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Estimate the pointwise derivatives at representative
# values of Ac
s <- avg_slopes(fit,
                variables = "Ac",
                newdata = datagrid(Ac = values,
                                   grid_type = "counterfactual"),
                by = "Ac")

# Plot the AMEF
ggplot(s, aes(x = Ac)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Ac", y = "dE[Y|A]/dA") +
  theme_bw()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data("msmdata")

Wmsm <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                             A_2 ~ X1_1 + X2_1 +
                               A_1 + X1_0 + X2_0,
                             A_3 ~ X1_2 + X2_2 +
                               A_2 + X1_1 + X2_1 +
                               A_1 + X1_0 + X2_0),
                        data = msmdata, method = "glm",
                        stabilize = TRUE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit <- glm_weightit(Y_B ~ A_1 * A_2 * A_3 * (X1_0 + X2_0),
                    data = msmdata, weightit = Wmsm,
                    family = binomial)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
(p <- avg_predictions(fit,
                      variables = c("A_1", "A_2", "A_3")))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
hypotheses(p, "reference")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Wm <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d,
               method = "cbps", estimand = "ATE",
               by = ~X5)

Wm

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit <- lm_weightit(Y_C ~ A * X5 * (X1 + X2 + X3),
                    data = d, weightit = Wm)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_comparisons(fit, variables = "A",
                by = "X5")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_comparisons(fit, variables = "A",
                by = "X5",
                hypothesis = "pairwise")

## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  avg_comparisons(fit, variables = "A",
#                  by = "X5",
#                  hypothesis = "reference") |>
#    hypotheses(joint = TRUE)

## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Generating data similar to Austin (2009) for demonstrating treatment effect estimation
#  gen_X <- function(n) {
#    X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
#    X[,5] <- as.numeric(X[,5] < .5)
#    X
#  }
#  
#  gen_Ac <- function(X) {
#    LP_A <- -1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
#    LP_A + rlogis(nrow(X))
#  }
#  
#  #~20% treated
#  gen_A <- function(Ac) {
#    1 * (Ac > 0)
#  }
#  
#  gen_Am <- function(A) {
#    factor(ifelse(A == 1, "T", sample(c("C1", "C2"), length(A), TRUE)))
#  }
#  
#  # Continuous outcome
#  gen_Y_C <- function(A, X) {
#    2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
#  }
#  #Conditional:
#  #  MD: 2
#  #Marginal:
#  #  MD: 2
#  
#  # Binary outcome
#  gen_Y_B <- function(A, X) {
#    LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
#    P_B <- plogis(LP_B)
#    rbinom(length(A), 1, P_B)
#  }
#  #Conditional:
#  #  OR:   2.4
#  #  logOR: .875
#  #Marginal:
#  #  RD:    .144
#  #  RR:   1.54
#  #  logRR: .433
#  #  OR:   1.92
#  #  logOR  .655
#  
#  # Survival outcome
#  gen_Y_S <- function(A, X) {
#    LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
#    sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
#  }
#  #Conditional:
#  #  HR:   2.4
#  #  logHR: .875
#  #Marginal:
#  #  HR:   1.57
#  #  logHR: .452
#  
#  set.seed(19599)
#  
#  n <- 2000
#  X <- gen_X(n)
#  Ac <- gen_Ac(X)
#  A <- gen_A(Ac)
#  Am <- gen_Am(A)
#  
#  Y_C <- gen_Y_C(A, X)
#  Y_B <- gen_Y_B(A, X)
#  Y_S <- gen_Y_S(A, X)
#  
#  d <- data.frame(A, Am, Ac, X, Y_C, Y_B, Y_S)

