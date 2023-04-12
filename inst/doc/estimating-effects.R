## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=T)
options(width = 200, digits= 4)

me_ok <- requireNamespace("marginaleffects", quietly = TRUE) &&
  requireNamespace("sandwich", quietly = TRUE)
su_ok <- requireNamespace("survival", quietly = TRUE)
boot_ok <- requireNamespace("boot", quietly = TRUE)

#Generating data similar to Austin (2009) for demonstrating treatment effect estimation
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

#~20% treated
gen_A <- function(X) {
  LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
  P_A <- plogis(LP_A)
  rbinom(nrow(X), 1, P_A)
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
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
A <- gen_A(X)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, X, Y_C, Y_B, Y_S)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(d)

## ----message=FALSE,warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("WeightIt")

## ----message=FALSE,warning=FALSE, eval=me_ok----------------------------------------------------------------------------------------------------------------------------------------------------------
library("marginaleffects")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#PS weighting for the ATE with a logistic regression PS
W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "glm", estimand = "ATE")
W

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Bring weights into the dataset
d$weights <- W$weights

#Linear model with covariates
fit <- lm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
           data = d, weights = weights)

## ---- eval=me_ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_comparisons(fit,
                variables = "A",
                vcov = "HC3",
                wts = "weights")

## ---- eval=me_ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
avg_predictions(fit,
                variables = "A",
                vcov = "HC3",
                wts = "weights")

## ---- eval=me_ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Logistic regression model with covariates
fit <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
           data = d, weights = weights,
           family = quasibinomial)

#Compute effects; RR and confidence interval
avg_comparisons(fit,
                variables = "A",
                vcov = "HC3",
                wts = "weights",
                comparison = "lnratioavg",
                transform = "exp")

## ---- eval=su_ok, message=F---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("survival")

#Cox Regression for marginal HR
coxph(Surv(Y_S) ~ A, data = d, weights = weights)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
d$Am <- factor(ifelse(d$A == 1, "T", sample(c("C1", "C2"), nrow(d), TRUE)))

table(d$Am)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "ebal", estimand = "ATT",
              focal = "T")
W

## ---- eval=me_ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Bring weights into the dataset
d$weights <- W$weights

#Fit the outcome model
fit <- lm(Y_C ~ Am * (X1 + X2 + X3 + X4 + X5 + 
                        X6 + X7 + X8 + X9),
          data = d, weights = weights)

#G-computation
p <- avg_predictions(fit,
                     variables = "Am",
                     vcov = "HC3",
                     newdata = subset(boot_data, A == 1),
                     wts = "weights")
p

hypotheses(p, "revpairwise")

## ---- eval = requireNamespace("survey", quietly = TRUE), message=F------------------------------------------------------------------------------------------------------------------------------------
library("survey")

#Declare a survey design using the estimated weights
des <- svydesign(~1, weights = ~weights, data = d)

#Fit the outcome model
fit <- svyglm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                           X6 + X7 + X8 + X9),
              design = des)

#G-computation for the difference in means
avg_comparisons(fit,
                variables = "A",
                wts = "(weights)")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #PS weighting for the ATE
  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9,
               data = boot_data,
               method = "glm", estimand = "ATT")
  
  #Bring weights into the dataset
  boot_data$weights <- W$weights
  
  #Fit outcome model
  fit <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9),
             data = boot_data, weights = weights,
             family = quasibinomial)
  
  #G-computation
  comp <- avg_comparisons(fit,
                          variables = "A",
                          vcov = FALSE,
                          newdata = subset(boot_data, A == 1),
                          wts = "weights",
                          comparison = "lnratioavg",
                          transform = "exp")
  
  comp$estimate
}

## ---- eval = boot_ok, message=F, warning=F------------------------------------------------------------------------------------------------------------------------------------------------------------
library("boot")
set.seed(54321)
boot_out <- boot(d, boot_fun, R = 199)

boot_out
boot.ci(boot_out, type = "perc")

## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
b <- {
  if (boot_ok) boot.ci(boot_out, type = "perc")
  else list(t0 = 1.522, percent = c(0, 0, 0, 1.305, 1.800))
}

## ---- eval=su_ok--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
boot_fun <- function(data, i) {
  boot_data <- data[i,]
  
  #PS weighting for the ATE
  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9,
               data = boot_data,
               method = "glm", estimand = "ATT")
  
  #Bring weights into the dataset
  boot_data$weights <- W$weights
  
  #Fit outcome model
  fit <- coxph(Surv(Y_S) ~ A, data = boot_data,
               weights = weights)
  
  #Return the coefficient on treatment
  coef(fit)["A"]
}

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Generating data similar to Austin (2009) for demonstrating treatment effect estimation
#  gen_X <- function(n) {
#    X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
#    X[,5] <- as.numeric(X[,5] < .5)
#    X
#  }
#  
#  #~20% treated
#  gen_A <- function(X) {
#    LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
#    P_A <- plogis(LP_A)
#    rbinom(nrow(X), 1, P_A)
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
#  A <- gen_A(X)
#  
#  Y_C <- gen_Y_C(A, X)
#  Y_B <- gen_Y_B(A, X)
#  Y_S <- gen_Y_S(A, X)
#  
#  d <- data.frame(A, X, Y_C, Y_B, Y_S)

