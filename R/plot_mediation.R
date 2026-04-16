#' Dose-Response Curve: Natural Direct and Indirect Effects
#'
#' @description
#' Plots NDE, NIE, and (optionally) total effect as smooth spline curves over
#' the full range of treatment values, with pointwise bootstrap confidence bands.
#' This is the signature visualisation of **RobustMediate** and is
#' publication-ready out of the box.
#'
#' @param x          A `robmedfit` object.
#' @param estimands  Character vector of estimands to display.
#'                   Any subset of `c("NDE", "NIE", "TE")`. Default `c("NDE","NIE")`.
#' @param show_total Shorthand for adding `"TE"` to `estimands`. Default `FALSE`.
#' @param facet      Logical. Split estimands into separate facets? Default `FALSE`.
#' @param ...        Ignored.
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, R = 50
#' )
#' plot_mediation(fit)
#' plot_mediation(fit, estimands = c("NDE", "NIE", "TE"), facet = TRUE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_ribbon geom_line
#'   scale_color_manual scale_fill_manual labs theme_minimal facet_wrap
#' @export
plot_mediation <- function(x, estimands = c("NDE","NIE"),
                           show_total = FALSE, facet = FALSE, ...) {

  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")

  if (show_total) estimands <- union(estimands, "TE")
  valid <- c("NDE","NIE","TE")
  bad   <- setdiff(estimands, valid)
  if (length(bad)) rlang::abort(paste0("Unknown estimand(s): ", paste(bad, collapse = ", ")))

  eff_df <- x$effects$curves
  eff_df <- eff_df[eff_df$estimand %in% estimands, , drop = FALSE]

  pal_colour <- c(NDE = "#534AB7", NIE = "#0F6E56", TE = "#993C1D")
  pal_fill   <- pal_colour

  p <- ggplot2::ggplot(
    eff_df,
    ggplot2::aes(x = dose, y = estimate,
                 colour = estimand, fill = estimand)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = x$meta$ref_dose, linetype = "dotted",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = pal_colour[estimands], name = NULL) +
    ggplot2::scale_fill_manual(  values = pal_fill[estimands],   name = NULL) +
    ggplot2::labs(
      x       = paste0("Treatment dose (", x$meta$treat_var, ")"),
      y       = "Effect estimate",
      caption = paste0("Dotted vertical line = reference dose (",
                       round(x$meta$ref_dose, 3),
                       ").  Shaded bands = ",
                       scales::percent(1 - x$meta$alpha), " bootstrap CI.")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (facet) p <- p + ggplot2::facet_wrap(~ estimand, scales = "free_y")

  p
}
