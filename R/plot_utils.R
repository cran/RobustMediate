# Suppress R CMD check NOTEs for ggplot2 aes() column-name variables
utils::globalVariables(c(
  # plot_mediation / compare_fits
  "dose", "estimate", "estimand", "lower", "upper", "model",

  # plot_balance
  "smd", "covariate", "timing",

  # plot_meditcv
  "corridor_lo", "corridor_hi", "r_crit", "r_obs",
  "x_pos", "exceeds_meditcv", "label_txt",

  # plot_meditcv_profile
  "delta", "r_attenuated", "pathway", "r_att",

  # plot_sensitivity
  "evalue", "rho", "effect",

  # plot_curvature
  "y",
  "sensitivity_meditcv", "plot_meditcv", "sensitivity_meditcv_profile"
))


#' Plot a robmedfit object
#'
#' @description
#' Dispatches to the appropriate plot function based on `type`.
#'
#' @param x    A `robmedfit` object.
#' @param type Which plot to produce. One of:
#'   \describe{
#'     \item{`"mediation"`}{Dose-response curve of NDE/NIE/TE (default).}
#'     \item{`"balance"`}{Dual love plot of covariate balance.}
#'     \item{`"sensitivity"`}{E-value x rho sensitivity contour.}
#'     \item{`"meditcv"`}{Pathway-specific medITCV robustness corridor.}
#'     \item{`"meditcv_profile"`}{medITCV robustness profile (fragility decomposition).}
#'     \item{`"curvature"`}{Dose-varying fragility (curvature-based sensitivity).}
#'   }
#' @param ...  Passed to the underlying plot function.
#' @return A `ggplot2` object.
#' @export
plot.robmedfit <- function(
    x,
    type = c("mediation", "balance", "sensitivity",
             "meditcv", "meditcv_profile", "curvature"),
    ...
) {
  type <- match.arg(type)

  switch(type,

         mediation   = plot_mediation(x, ...),
         balance     = plot_balance(x, ...),
         sensitivity = plot_sensitivity(x, ...),

         meditcv = {
           if (is.null(x$meditcv)) {
             rlang::abort("medITCV not found in fit object. Re-run robustmediate().")
           }
           plot_meditcv(x$meditcv, ...)
         },

         meditcv_profile = {
           obj <- x$meditcv_profile
           if (is.null(obj)) {
             rlang::abort("medITCV robustness profile not found in fit. Re-run robustmediate().")
           }
           plot_meditcv_profile(obj, ...)
         },

         curvature = {
           curv <- sensitivity_curvature(x, ...)
           plot_curvature(curv, ref_dose = x$meta$ref_dose)
         }
  )
}


#' Compare Two robmedfit Objects Side by Side
#'
#' @description
#' Overlays the NDE/NIE/TE curves from two `robmedfit` objects on the same
#' panel. Useful for sensitivity comparisons (e.g. different spline degrees,
#' trimming thresholds, or model specifications).
#'
#' @param fit1      First `robmedfit` object.
#' @param fit2      Second `robmedfit` object.
#' @param label1    Label for `fit1`. Default `"Model 1"`.
#' @param label2    Label for `fit2`. Default `"Model 2"`.
#' @param estimands Estimands to display. Default `c("NDE","NIE")`.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' fit_a <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, spline_df = 3, R = 50
#' )
#' fit_b <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, spline_df = 6, R = 50
#' )
#' compare_fits(fit_a, fit_b, label1 = "df=3", label2 = "df=6")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_hline
#'   scale_color_manual scale_fill_manual scale_linetype_manual labs
#'   theme_minimal theme
#' @export
compare_fits <- function(fit1, fit2,
                         label1    = "Model 1",
                         label2    = "Model 2",
                         estimands = c("NDE", "NIE")) {
  if (!inherits(fit1, "robmedfit") || !inherits(fit2, "robmedfit")) {
    rlang::abort("Both arguments must be robmedfit objects.")
  }

  mk_df <- function(fit, lbl) {
    d       <- fit$effects$curves
    d       <- d[d$estimand %in% estimands, , drop = FALSE]
    d$model <- lbl
    d
  }
  combined       <- rbind(mk_df(fit1, label1), mk_df(fit2, label2))
  combined$model <- factor(combined$model, levels = c(label1, label2))

  pal <- c(NDE = "#534AB7", NIE = "#0F6E56", TE = "#993C1D")

  ggplot2::ggplot(combined,
                  ggplot2::aes(x = dose, y = estimate,
                               colour = estimand, fill = estimand, linetype = model)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.10, colour = NA) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::scale_colour_manual(values = pal[estimands], name = "Estimand") +
    ggplot2::scale_fill_manual(  values = pal[estimands], name = "Estimand") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = "Model") +
    ggplot2::labs(
      x       = paste0("Treatment dose (", fit1$meta$treat_var, ")"),
      y       = "Effect estimate",
      caption = paste0("Solid = ", label1, "  |  Dashed = ", label2,
                       ".  Shaded bands = bootstrap CI.")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")
}


#' Extract Tipping-Point Table
#'
#' @description
#' Returns a formatted data frame of sensitivity tipping points: the minimum
#' E-value and minimum |rho| required to nullify the NIE. Designed for
#' direct insertion into a table in a manuscript.
#'
#' @param x A `robmedfit` object.
#'
#' @return A data frame with columns `parameter`, `tipping_value`,
#'   `interpretation`.
#'
#' @export
tipping_table <- function(x) {
  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")
  tip <- x$sensitivity$tipping
  eff <- x$effects$summary
  data.frame(
    parameter      = c("E-value (NIE)", "Min |rho| (NIE)"),
    tipping_value  = c(round(tip$evalue_NIE, 2), round(abs(tip$rho_min), 2)),
    interpretation = c(
      paste0("Unmeasured confounder must be associated ",
             round(tip$evalue_NIE, 1), "x with both treatment & outcome"),
      paste0("|rho| must exceed ", round(abs(tip$rho_min), 2),
             " to nullify NIE = ", round(eff$NIE, 3))
    ),
    stringsAsFactors = FALSE
  )
}
