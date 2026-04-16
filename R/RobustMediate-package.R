#' RobustMediate: Causal Mediation Analysis with Diagnostics and Sensitivity Analysis
#'
#' @description
#' **RobustMediate** provides a workflow for causal mediation analysis with
#' continuous treatments using inverse probability weighting (IPW), diagnostic
#' tools, and sensitivity analysis.
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
#' - **`plot_sensitivity()`** - Displays a bivariate sensitivity surface based on
#'   E-values and sequential ignorability violations parameterized by `rho`.
#'
#' - **`sensitivity_meditcv()`** - Computes pathway-specific mediation ITCV
#'   (`medITCV`) diagnostics based on Frank's impact threshold for a
#'   confounding variable framework.
#'
#' - **`plot_meditcv()`** - Displays pathway-specific medITCV robustness
#'   corridors for the a-path and b-path.
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
#' plot(fit, type = "meditcv")
#' diagnose(fit)
#' ```
#'
#' @section Sensitivity interpretation:
#' The E-value x `rho` surface is a bivariate robustness display rather than a
#' single unified causal model. It is intended to help users examine how large
#' different classes of unmeasured confounding would need to be to attenuate or
#' nullify the estimated indirect effect.
#'
#' The mediation ITCV (`medITCV`) is reported separately for the a-path and
#' b-path. The indirect-effect summary is interpreted as a minimum-path
#' robustness bound governed by the weaker pathway.
#'
#' @references
#' Frank, K. A. (2000). Impact of a confounding variable on a regression
#' coefficient. *Sociological Methods & Research*, 29(2), 147--194.
#'
#' Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference, and
#' sensitivity analysis for causal mediation effects. *Psychological Methods*,
#' 15(4), 309--334.
#'
#' VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in observational
#' research: Introducing the E-value. *Annals of Internal Medicine*, 167(4),
#' 268--274.
#'
#' @docType package
#' @name RobustMediate-package
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "RobustMediate ", utils::packageVersion("RobustMediate"),
    " - Causal mediation analysis with diagnostics and sensitivity analysis.\n",
    "  Docs: ?robustmediate | ?plot_balance | ?plot_mediation | ?plot_sensitivity | ?plot_meditcv"
  )
}
