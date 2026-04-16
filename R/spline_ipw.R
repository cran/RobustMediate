#' Spline-Based Generalised Propensity Score for Continuous Treatments
#'
#' @description
#' Internal machinery that fits a **natural spline** treatment model to estimate
#' the generalised propensity score (GPS) for a continuous exposure. This
#' produces more flexible, positivity-respecting stabilised IPW weights than a
#' simple linear model.
#'
#' Users do not call these functions directly — they are invoked automatically
#' by `robustmediate()` when `spline_df > 1`.
#'
#' @name spline_ipw
#' @keywords internal
NULL


## ── Spline GPS model ─────────────────────────────────────────────────────────

#' Fit a natural-spline treatment model and return the GPS model object
#'
#' @param treatment_formula Original treatment formula (e.g., `X ~ Z1 + Z2`).
#' @param data              Data frame.
#' @param df                Degrees of freedom for [splines::ns()]. Default 4.
#' @param family            GLM family. Default `gaussian()`.
#'
#' @return A fitted `glm` object with a spline basis for the intercept term
#'   replaced by `ns(X, df)` on the *response* side — i.e. the model
#'   `X_spline ~ Z` where `X_spline` is the ns-expanded treatment.
#'   In practice we use the approach: model X | Z with a Gaussian GLM using
#'   the raw X, then use the **residual SD** from that model to form GPS
#'   densities. The spline enters on the *covariate* side of mediator/outcome
#'   models for the dose-response curve.
#'
#' @keywords internal
.fit_gps_model <- function(treatment_formula, data, df = 4,
                            family = stats::gaussian()) {
  # Augment RHS covariates with spline basis for any numeric predictors
  # that appear non-linearly related to treatment.
  # Strategy: fit a standard GLM; check residual patterns; if df > 1, refit
  # with ns() on the treatment variable itself isn't the right move —
  # instead we use ns() on continuous *covariates* to capture nonlinear
  # confounding.

  rhs_vars  <- all.vars(treatment_formula)[-1]
  treat_var <- all.vars(treatment_formula)[1]

  # Identify continuous RHS variables (numeric, > 10 unique values)
  continuous_rhs <- Filter(function(v) {
    is.numeric(data[[v]]) && length(unique(data[[v]])) > 10
  }, rhs_vars)

  if (df > 1 && length(continuous_rhs) > 0) {
    # Build augmented formula: replace each continuous RHS variable with ns()
    new_terms <- lapply(rhs_vars, function(v) {
      if (v %in% continuous_rhs) {
        paste0("splines::ns(", v, ", df = ", df, ")")
      } else {
        v
      }
    })
    new_formula <- stats::as.formula(
      paste(treat_var, "~", paste(new_terms, collapse = " + "))
    )
  } else {
    new_formula <- treatment_formula
  }

  stats::glm(new_formula, data = data, family = family)
}


#' Compute stabilised IPW weights from a GPS model (continuous treatment)
#'
#' @description
#' Implements the stabilised weight formula from Robins, Hernán & Brumback
#' (2000) for continuous exposures:
#'
#' \deqn{w_i = \frac{f(X_i \mid \bar{X}, \hat{\sigma}^2_m)}{f(X_i \mid Z_i, \hat{\sigma}^2_c)}}
#'
#' where the numerator is the marginal density of X and the denominator is
#' the conditional density given covariates. Both are evaluated as Gaussian
#' PDFs using the model residual standard deviation.
#'
#' Weights are trimmed at the 1st and 99th percentile to reduce the influence
#' of extreme propensity scores (positivity violations).
#'
#' @param gps_model  Fitted GPS model from `.fit_gps_model()`.
#' @param data       Original data frame.
#' @param treat_var  Name of the treatment variable (character).
#' @param trim       Quantile trimming bounds. Default `c(0.01, 0.99)`.
#'
#' @return Numeric vector of stabilised, trimmed IPW weights (length = nrow(data)).
#'
#' @references
#' Robins, J. M., Hernán, M. A., & Brumback, B. (2000). Marginal structural
#' models and causal inference in epidemiology. *Epidemiology*, 11(5), 550–560.
#'
#' Hirano, K. & Imbens, G. W. (2004). The propensity score with continuous
#' treatments. In A. Gelman & X.-L. Meng (Eds.), *Applied Bayesian modeling
#' and causal inference from incomplete-data perspectives* (pp. 73–84). Wiley.
#'
#' @keywords internal
.compute_gps_weights <- function(gps_model, data, treat_var, trim = c(0.01, 0.99)) {
  X_obs  <- data[[treat_var]]
  mu_hat <- stats::fitted(gps_model)
  # Conditional SD from model residuals
  sigma_c <- stats::sd(stats::residuals(gps_model))

  # Marginal density (numerator): N(X | mean(X), sd(X))
  mu_m   <- mean(X_obs)
  sigma_m <- stats::sd(X_obs)
  f_num  <- stats::dnorm(X_obs, mean = mu_m,   sd = sigma_m)

  # Conditional density (denominator): N(X | mu_hat(Z), sigma_c)
  f_den  <- stats::dnorm(X_obs, mean = mu_hat, sd = sigma_c)
  f_den  <- pmax(f_den, 1e-8)   # floor to avoid division by zero

  w <- f_num / f_den

  # Percentile trimming
  lo <- stats::quantile(w, trim[1])
  hi <- stats::quantile(w, trim[2])
  w  <- pmin(pmax(w, lo), hi)

  w
}


#' Effective Sample Size for IPW weights
#'
#' @description
#' Returns the effective sample size (ESS) of the weighted sample:
#' \eqn{ESS = (\sum w_i)^2 / \sum w_i^2}
#' Values well below n indicate severe weighting and likely positivity problems.
#'
#' @param weights Numeric vector of IPW weights.
#' @return Scalar ESS value.
#'
#' @keywords internal
.ess <- function(weights) {
  sum(weights)^2 / sum(weights^2)
}


#' Summarise GPS weight distribution
#'
#' @description
#' Prints a compact diagnostic table for IPW weights: range, mean, ESS, and
#' a flag if ESS < 0.4 * n (severe imbalance warning).
#'
#' @param weights   Numeric weight vector.
#' @param n         Original sample size.
#' @param pathway   Label for the pathway (e.g. `"treatment"`).
#'
#' @return Invisibly returns a named list with `min`, `max`, `mean`, `ess`, `ess_ratio`.
#'
#' @keywords internal
.summarise_weights <- function(weights, n, pathway = "treatment") {
  ess       <- .ess(weights)
  ess_ratio <- ess / n

  if (ess_ratio < 0.4) {
    rlang::warn(paste0(
      "GPS weights for the ", pathway, " pathway have ESS ratio = ",
      round(ess_ratio, 2), " (< 0.40). ",
      "This suggests possible positivity violations. ",
      "Consider increasing spline_df, trimming more aggressively, ",
      "or inspecting overlap with plot_balance()."
    ))
  }

  invisible(list(
    min       = min(weights),
    max       = max(weights),
    mean      = mean(weights),
    ess       = ess,
    ess_ratio = ess_ratio
  ))
}


## ── Dose-grid prediction helpers ─────────────────────────────────────────────

#' Build a prediction data frame for a specific dose value
#'
#' @description
#' Returns a copy of `data` with the treatment variable set to `dose`, used
#' to predict counterfactual mediator and outcome means under intervention.
#'
#' @param data      Original data frame.
#' @param treat_var Treatment variable name.
#' @param dose      Scalar dose value.
#'
#' @return Data frame with `treat_var` replaced by `dose`.
#'
#' @keywords internal
.dose_df <- function(data, treat_var, dose) {
  d <- data
  d[[treat_var]] <- dose
  d
}


#' Predict mean potential outcome under a dose pair (x, x_ref)
#'
#' @description
#' Implements the g-computation identification formula for the NDE and NIE:
#'
#' \deqn{NDE(x, x^*) = E[Y(x, M(x^*))] - E[Y(x^*, M(x^*))]}
#' \deqn{NIE(x, x^*) = E[Y(x, M(x))]   - E[Y(x,  M(x^*))]}
#'
#' The mediator is predicted at the **individual** level (not as a scalar mean),
#' correctly marginalising over each person's covariate values.
#'
#' @param models    List with elements `treatment`, `mediator`, `outcome`.
#' @param data      Original data frame.
#' @param x         Dose value (the "treated" level).
#' @param x_ref     Reference dose (counterfactual baseline).
#' @param treat_var Character. Name of treatment variable.
#' @param med_var   Character. Name of mediator variable.
#' @param out_var   Character. Name of outcome variable.
#'
#' @return Named numeric vector: `NDE`, `NIE`, `TE`.
#'
#' @keywords internal
.predict_effects <- function(models, data, x, x_ref,
                              treat_var, med_var, out_var) {
  d_x   <- .dose_df(data, treat_var, x)
  d_ref <- .dose_df(data, treat_var, x_ref)

  # Individual-level predicted mediators (vectors, not scalars)
  M_hat_x   <- stats::predict(models$mediator, newdata = d_x,   type = "response")
  M_hat_ref <- stats::predict(models$mediator, newdata = d_ref, type = "response")

  # Three cross-world data frames
  d_x_Mx     <- d_x;   d_x_Mx[[med_var]]    <- M_hat_x    # X=x,  M=M(x)
  d_ref_Mref <- d_ref; d_ref_Mref[[med_var]] <- M_hat_ref  # X=x*, M=M(x*)
  d_x_Mref   <- d_x;   d_x_Mref[[med_var]]   <- M_hat_ref  # X=x,  M=M(x*)

  Y_x_Mx     <- mean(stats::predict(models$outcome, newdata = d_x_Mx,     type = "response"))
  Y_ref_Mref <- mean(stats::predict(models$outcome, newdata = d_ref_Mref, type = "response"))
  Y_x_Mref   <- mean(stats::predict(models$outcome, newdata = d_x_Mref,   type = "response"))

  NDE <- Y_x_Mref  - Y_ref_Mref   # E[Y(x, M(x*))] - E[Y(x*, M(x*))]
  NIE <- Y_x_Mx    - Y_x_Mref     # E[Y(x, M(x))]  - E[Y(x, M(x*))]
  TE  <- Y_x_Mx    - Y_ref_Mref   # = NDE + NIE

  c(NDE = NDE, NIE = NIE, TE = TE)
}
