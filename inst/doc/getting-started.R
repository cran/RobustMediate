## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.5
)

## ----data---------------------------------------------------------------------
set.seed(42)
n  <- 600
Z1 <- rnorm(n)
Z2 <- rbinom(n, 1, 0.5)
X  <- 0.5 * Z1 + 0.3 * Z2 + rnorm(n)  # continuous treatment
M  <- 0.4 * X  + 0.2 * Z1 + rnorm(n)  # mediator
Y  <- 0.3 * X  + 0.5 * M  + 0.1 * Z1 + rnorm(n)  # outcome
dat <- data.frame(Y, X, M, Z1, Z2)

## ----fit, message = FALSE-----------------------------------------------------
library(RobustMediate)

fit <- robustmediate(
  treatment_formula = X ~ Z1 + Z2,
  mediator_formula  = M ~ X + Z1 + Z2,
  outcome_formula   = Y ~ X + M + Z1 + Z2,
  data    = dat,
  R       = 50,    # use >= 500 in practice
  verbose = FALSE
)
print(fit)

## ----love-plot----------------------------------------------------------------
plot_balance(fit)

## ----mediation-plot-----------------------------------------------------------
plot_mediation(fit, estimands = c("NDE", "NIE", "TE"))

## ----sensitivity-plot---------------------------------------------------------
sensitivity_rho(fit)
plot_sensitivity(fit)

## ----diagnose-----------------------------------------------------------------
diagnose(fit)

## ----paths--------------------------------------------------------------------
sp <- sensitivity_paths(fit)
sp
plot_paths(sp)

## ----cluster, eval = FALSE----------------------------------------------------
# data(sim_mediation)
# fit_clus <- robustmediate(
#   treatment_formula = X ~ Z1 + Z2 + Z3,
#   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#   data        = sim_mediation,
#   cluster_var = "school",
#   R           = 500
# )

