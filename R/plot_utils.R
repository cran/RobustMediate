# Suppress R CMD check NOTEs for ggplot2 aes() column-name variables
utils::globalVariables(c(
  # plot_mediation / compare_fits
  "dose", "estimate", "estimand", "lower", "upper", "model",
  # plot_balance
  "smd", "covariate", "timing",
  # plot_sensitivity
  "rho", "nie", "stat", "line",
  # plot_curvature
  "y"
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
#'     \item{`"sensitivity"`}{Indirect effect as a function of `rho`
#'       (see [sensitivity_rho()]).}
#'     \item{`"curvature"`}{Effect curve and its curvature (see
#'       [effect_curvature()]).}
#'   }
#' @param ...  Passed to the underlying plot function.
#' @return A `ggplot2` object.
#' @export
plot.robmedfit <- function(
    x,
    type = c("mediation", "balance", "sensitivity", "curvature"),
    ...
) {
  type <- match.arg(type)

  switch(type,

         mediation   = plot_mediation(x, ...),
         balance     = plot_balance(x, ...),
         sensitivity = plot_sensitivity(x, ...),

         curvature = {
           curv <- effect_curvature(x, ...)
           plot_curvature(curv, ref_dose = x$meta$ref_dose)
         }
  )
}


#' Compare Two robmedfit Objects Side by Side
#'
#' @description
#' Overlays the NDE/NIE/TE curves from two `robmedfit` objects on the same
#' panel. Useful for comparing model specifications (e.g. different
#' covariate sets).
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
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' fit_b <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2,
#'   mediator_formula  = M ~ X + Z1 + Z2,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2,
#'   data = sim_mediation, R = 10, verbose = FALSE,
#'   dose_grid = seq(-2, 2, length.out = 25)
#' )
#' compare_fits(fit_a, fit_b, label1 = "With Z3", label2 = "Without Z3")
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
