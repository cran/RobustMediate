#' Print a robmedfit object
#'
#' @param x A `robmedfit` object returned by [robustmediate()].
#' @param ... Ignored.
#'
#' @return The input object, invisibly.
#' @export
print.robmedfit <- function(x, ...) {
  cat("-- RobustMediate fit ------------------------------------------\n")
  cat(sprintf(
    "  Treatment: %s  |  Mediator: %s  |  Outcome: %s\n",
    x$meta$treat_var, x$meta$med_var, x$meta$out_var
  ))
  cat(sprintf(
    "  N = %d  |  Ref dose = %.3f  |  Bootstrap reps = %d\n",
    x$meta$n_obs, x$meta$ref_dose, x$meta$R
  ))

  eff <- x$effects$summary
  cat(
    "\n  Effects at focal dose (",
    round(eff$focal_dose, 3),
    " vs ref ",
    round(x$meta$ref_dose, 3),
    "):\n",
    sep = ""
  )

  for (est in c("NDE", "NIE", "TE")) {
    cat(sprintf(
      "    %-4s  %.4f  [%.4f, %.4f]\n",
      est,
      eff[[est]],
      eff[[paste0(est, "_lo")]],
      eff[[paste0(est, "_hi")]]
    ))
  }

  cat(sprintf("    %% mediated: %.1f%%\n", eff$pct_mediated))

  bal <- x$balance$summary_stats
  cat("\n  Balance (max |SMD| after weighting):\n")
  cat(sprintf(
    "    Treatment pathway: %.3f  (%d covariate(s) above 0.10)\n",
    bal$treatment$max_smd, bal$treatment$n_above
  ))
  cat(sprintf(
    "    Mediator pathway:  %.3f  (%d covariate(s) above 0.10)\n",
    bal$mediator$max_smd, bal$mediator$n_above
  ))

  if (!is.null(x$sensitivity)) {
    cat(sprintf(
      "\n  Sensitivity: indirect effect is zero at rho = %.3f (see sensitivity_rho())\n",
      x$sensitivity$rho_zero
    ))
  } else {
    cat("\n  Sensitivity: not available (requires linear mediator and outcome models)\n")
  }

  cat("---------------------------------------------------------------\n")
  invisible(x)
}


#' Summary method for robmedfit objects
#'
#' @param object A `robmedfit` object.
#' @param ... Ignored.
#'
#' @return A list with effect, balance, and sensitivity summaries, invisibly.
#' @export
summary.robmedfit <- function(object, ...) {
  print(object)

  out <- list(
    effects = object$effects$summary,
    balance = object$balance$summary_stats,
    sensitivity = object$sensitivity
  )

  invisible(out)
}
