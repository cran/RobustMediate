#' Love Plot: Balance Diagnostics for Both Pathways
#'
#' @description
#' Produces a publication-ready love plot showing standardised mean differences
#' (SMDs) before and after IPW weighting for **both** the treatment and mediator
#' pathways — stacked vertically in a single panel. This dual-pathway display is
#' unique to **RobustMediate**; no other mediation package provides it.
#'
#' @param x         A `robmedfit` object.
#' @param threshold Absolute SMD threshold displayed as dashed reference lines.
#'                  Reviewers conventionally accept |SMD| < 0.10. Default `0.1`.
#' @param pathways  Character vector indicating which pathways to show.
#'                  Options: `"treatment"`, `"mediator"`, or both (default).
#' @param ...       Ignored (for S3 consistency).
#'
#' @return A `ggplot2` object. Add layers or themes as usual.
#'
#' @examples
#' \donttest{
#' fit <- robustmediate(
#'   treatment_formula = X ~ Z1 + Z2 + Z3,
#'   mediator_formula  = M ~ X + Z1 + Z2 + Z3,
#'   outcome_formula   = Y ~ X + M + Z1 + Z2 + Z3,
#'   data = sim_mediation, R = 50
#' )
#' plot_balance(fit)
#' plot_balance(fit, threshold = 0.05, pathways = "treatment")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_point facet_wrap labs
#'   theme_minimal theme element_line
#' @importFrom stats aggregate
#' @export
plot_balance <- function(x, threshold = 0.1,
                         pathways = c("treatment", "mediator"), ...) {

  if (!inherits(x, "robmedfit")) rlang::abort("`x` must be a robmedfit object.")

  bal_df <- x$balance$summary
  bal_df <- bal_df[bal_df$pathway %in% pathways, , drop = FALSE]

  # Reshape to long for pre/post colour mapping
  long <- rbind(
    data.frame(bal_df[, c("covariate","pathway")], timing = "Before", smd = bal_df$smd_pre),
    data.frame(bal_df[, c("covariate","pathway")], timing = "After",  smd = bal_df$smd_post)
  )
  long$timing <- factor(long$timing, levels = c("Before","After"))

  # Order covariates by post-weighting |SMD| descending
  ord <- aggregate(abs(smd) ~ covariate, data = long[long$timing == "After", ], FUN = max)
  ord <- ord[order(-ord[["abs(smd)"]] ), "covariate"]
  long$covariate <- factor(long$covariate, levels = rev(ord))

  ggplot2::ggplot(long, ggplot2::aes(x = smd, y = covariate,
                                     colour = timing, shape = timing)) +
    ggplot2::geom_vline(xintercept = c(-threshold, threshold),
                        linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_colour_manual(
      values = c(Before = "#D85A30", After = "#1D9E75"),
      name   = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c(Before = 1, After = 19),
      name   = NULL
    ) +
    ggplot2::facet_wrap(~ pathway, ncol = 1,
                        labeller = ggplot2::labeller(
                          pathway = c(treatment = "Treatment pathway",
                                      mediator  = "Mediator pathway")
                        )) +
    ggplot2::labs(
      x       = "Standardised mean difference",
      y       = NULL,
      caption = paste0("Dashed lines at +/-", threshold,
                       ".  Open circles = before weighting; filled = after.")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.2),
      legend.position    = "top"
    )
}
