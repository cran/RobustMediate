# Re-export broom generics so users can call tidy(), glance(), augment()
# without needing to attach broom themselves.

#' @importFrom broom tidy
#' @export
broom::tidy

#' @importFrom broom glance
#' @export
broom::glance

#' @importFrom broom augment
#' @export
broom::augment


#' Tidy a robmedfit object (broom-compatible)
#'
#' @description
#' Returns a tidy data frame of the mediation effect estimates (NDE, NIE, TE)
#' at the reference dose, with confidence intervals. Compatible with
#' `broom::tidy()` and the broader `tidymodels` ecosystem.
#'
#' @param x       A `robmedfit` object.
#' @param conf.int Logical. Include confidence interval columns? Default `TRUE`.
#' @param ...     Ignored.
#'
#' @return A data frame with columns `term`, `estimate`, `conf.low`,
#'   `conf.high`, and `ref_dose`.
#'
#' @examples
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, R = 50
#' )
#' tidy(fit)
#' }
#'
#' @export
tidy.robmedfit <- function(x, conf.int = TRUE, ...) {
  eff <- x$effects$summary
  out <- data.frame(
    term     = c("NDE",    "NIE",    "TE"),
    estimate = c(eff$NDE,  eff$NIE,  eff$TE),
    ref_dose = x$meta$ref_dose,
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    out$conf.low  <- c(eff$NDE_lo, eff$NIE_lo, eff$TE_lo)
    out$conf.high <- c(eff$NDE_hi, eff$NIE_hi, eff$TE_hi)
  }
  out
}


#' Glance at a robmedfit object (broom-compatible)
#'
#' @description
#' Returns a one-row summary of the fit: sample size, bootstrap reps, reference
#' dose, percentage mediated, and the two tipping-point sensitivity values.
#'
#' @param x   A `robmedfit` object.
#' @param ... Ignored.
#'
#' @return A one-row data frame.
#'
#' @examples
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, R = 50
#' )
#' glance(fit)
#' }
#'
#' @export
glance.robmedfit <- function(x, ...) {
  data.frame(
    n_obs         = x$meta$n_obs,
    R             = x$meta$R,
    ref_dose      = x$meta$ref_dose,
    pct_mediated  = x$effects$summary$pct_mediated,
    evalue_NIE    = x$sensitivity$tipping$evalue_NIE,
    rho_min_NIE   = x$sensitivity$tipping$rho_min,
    max_smd_treat = x$balance$summary_stats$treatment$max_smd,
    max_smd_med   = x$balance$summary_stats$mediator$max_smd,
    stringsAsFactors = FALSE
  )
}


#' Augment data with fitted values from a robmedfit object (broom-compatible)
#'
#' @description
#' Returns the original data augmented with IPW weights, fitted mediator values,
#' and fitted outcome values from the pathway models.
#'
#' @param x    A `robmedfit` object.
#' @param data Optional data frame to augment. Defaults to the data stored
#'             in the outcome model's `$data` slot.
#' @param ...  Ignored.
#'
#' @return The original data frame with additional columns:
#' \describe{
#'   \item{`.ipw_weight`}{Stabilised inverse probability weights.}
#'   \item{`.fitted_mediator`}{Predicted mediator values.}
#'   \item{`.fitted_outcome`}{Predicted outcome values.}
#'   \item{`.resid_mediator`}{Residuals from the mediator model.}
#'   \item{`.resid_outcome`}{Residuals from the outcome model.}
#' }
#'
#' @examples
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, R = 50
#' )
#' aug <- augment(fit)
#' hist(aug$.ipw_weight)
#' }
#'
#' @export
augment.robmedfit <- function(x, data = NULL, ...) {
  if (is.null(data)) {
    data <- tryCatch(x$models$outcome$data, error = function(e) NULL)
    if (is.null(data)) {
      rlang::abort("Could not extract data from model. Please supply `data` explicitly.")
    }
  }
  ipw <- .compute_gps_weights(x$models$treatment, data, x$meta$treat_var)

  data$.ipw_weight      <- ipw
  data$.fitted_mediator <- stats::fitted(x$models$mediator)
  data$.fitted_outcome  <- stats::fitted(x$models$outcome)
  data$.resid_mediator  <- stats::residuals(x$models$mediator)
  data$.resid_outcome   <- stats::residuals(x$models$outcome)
  data
}


#' Coerce effects curve to a data frame
#'
#' @description
#' Convenience accessor returning the NDE/NIE/TE curve data frame.
#'
#' @param x   A `robmedfit` object.
#' @param ... Ignored.
#'
#' @return A data frame with columns `dose`, `estimand`, `estimate`,
#'   `lower`, `upper`.
#'
#' @export
as.data.frame.robmedfit <- function(x, ...) {
  x$effects$curves
}
