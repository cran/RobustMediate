#' RobustMediate: Causal Mediation Analysis with Diagnostics and Sensitivity Analysis
#'
#' @description
#' **RobustMediate** provides a workflow for causal mediation analysis with
#' continuous treatments: g-computation estimates of natural direct and
#' indirect effects, balance diagnostics, and sensitivity analysis.
#'
#' Main functions include:
#'
#' - **`robustmediate()`** - Fits treatment, mediator, and outcome models and
#'   stores precomputed results for downstream plotting and reporting.
#'
#' - **`plot_balance()`** - Displays covariate balance before and after weighting
#'   for both the treatment and mediator pathways using standardized mean
#'   differences.
#'
#' - **`plot_mediation()`** - Plots estimated natural direct effects (NDE) and
#'   natural indirect effects (NIE) over the treatment range, with pointwise
#'   uncertainty bands.
#'
#' - **`sensitivity_rho()`** and **`plot_sensitivity()`** - Sensitivity of the
#'   indirect effect to confounding of the mediator-outcome relationship, using
#'   the correlation parameter `rho` of Imai, Keele, and Yamamoto (2010).
#'
#' - **`diagnose()`** - Produces a formatted diagnostic summary of balance,
#'   mediation effects, and sensitivity results.
#'
#' @section Getting started:
#' ```r
#' library(RobustMediate)
#'
#' data(sim_mediation)
#'
#' fit <- robustmediate(
#'   X ~ Z1 + Z2 + Z3,
#'   M ~ X + Z1 + Z2 + Z3,
#'   Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation,
#'   R = 500
#' )
#'
#' plot(fit)
#' plot(fit, type = "balance")
#' plot(fit, type = "sensitivity")
#' diagnose(fit)
#' ```
#'
#' @section Sensitivity interpretation:
#' `sensitivity_rho()` reports how the indirect effect changes as the
#' correlation `rho` between the mediator-model and outcome-model errors moves
#' away from zero, and the value of `rho` at which it reaches zero. It applies
#' to linear models only and assumes the treatment is ignorable given the
#' covariates.
#'
#' @references
#' Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference and
#' sensitivity analysis for causal mediation effects. *Statistical Science*,
#' 25(1), 51--71.
#'
#' @docType package
#' @name RobustMediate-package
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "RobustMediate ", utils::packageVersion("RobustMediate"),
    " - Causal mediation analysis with diagnostics and sensitivity analysis.\n",
    "  Docs: ?robustmediate | ?plot_balance | ?plot_mediation | ?sensitivity_rho"
  )
}
