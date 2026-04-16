#' Sensitivity Contour Plot: Joint E-value x Rho Display
#'
#' @description
#' Renders the novel bivariate robustness map unique to **RobustMediate**: a 2-D
#' heatmap where the x-axis is the E-value (VanderWeele-style unmeasured
#' treatment–outcome confounding) and the y-axis is Imai's sequential-ignorability
#' violation parameter rho. Contour lines show where the mediation effect crosses
#' zero, so readers can judge robustness to **two different** sensitivity
#' dimensions simultaneously.
#'
#' This visualisation does not exist elsewhere in the R ecosystem. The correct
#' interpretation is as a *bivariate robustness display*, not a joint causal
#' model — see the package paper for theoretical justification.
#'
#' @param x             A `robmedfit` object.
#' @param annotate_zero Logical. Draw a bold dashed zero-crossing contour?
#'                      Default `TRUE`. Highly recommended for applied papers.
#' @param n_breaks      Number of fill colour breaks in the heatmap. Default `12`.
#' @param palette       RColorBrewer palette name. Default `"RdYlGn"`.
#' @param ...           Ignored.
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
#' plot_sensitivity(fit)
#' plot_sensitivity(fit, annotate_zero = FALSE, palette = "PuOr")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_contour_filled geom_contour
#'   scale_fill_brewer labs theme_minimal
#' @export
plot_sensitivity <- function(x, annotate_zero = TRUE,
                             n_breaks = 12, palette = "RdYlGn", ...) {

  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")

  sens_df <- x$sensitivity$surface
  tip     <- x$sensitivity$tipping

  brks <- seq(min(sens_df$effect, na.rm = TRUE),
              max(sens_df$effect, na.rm = TRUE),
              length.out = n_breaks)

  p <- ggplot2::ggplot(sens_df, ggplot2::aes(x = evalue, y = rho, z = effect)) +
    ggplot2::geom_contour_filled(breaks = brks) +
    ggplot2::scale_fill_brewer(palette = palette, direction = 1,
                               name = "Mediation\neffect (NIE)") +
    ggplot2::labs(
      x       = "E-value  (unmeasured treatment\u2013outcome confounding)",
      y       = expression(rho ~ "(sequential ignorability violation)"),
      caption = paste0(
        "Zero-crossing: E-value ~= ", round(tip$evalue_NIE, 2),
        " | min |rho| ~= ", round(abs(tip$rho_min), 2),
        ".\nInterpreted as a bivariate robustness display (not a joint causal model)."
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right")

  if (annotate_zero) {
    p <- p + ggplot2::geom_contour(
      breaks    = 0,
      colour    = "black",
      linewidth = 0.8,
      linetype  = "dashed"
    )
  }

  p
}
