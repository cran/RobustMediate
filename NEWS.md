# RobustMediate 0.3.0

## New features

* `sensitivity_paths()` reports robustness values for the a-path, the
  b-path, and the indirect effect in the partial R-squared metric of Cinelli
  and Hazlett (2020). One confounder enters both paths through its partial
  R-squared with the mediator. For the indirect effect it gives the robustness
  value under the joint-significance test (equal to the smaller path value)
  and under the Sobel test (always smaller), and their ratio. It accepts a
  `robmedfit` object or fitted `lm` mediator and outcome models.

* `adjust_paths()` returns the adjusted path and Sobel statistics for a
  confounder with separate partial R-squared values for the treatment,
  mediator, and outcome.

* `plot_paths()` plots the adjusted statistics against confounder strength.

* `effect_curvature()` replaces `sensitivity_curvature()` (now deprecated).
  It describes the shape of an effect curve and is not a sensitivity analysis.
  It returns curvature exactly zero when the curve is linear by construction
  (previously it plotted rounding error), and gives percentile bootstrap
  intervals for the curvature otherwise. `plot_curvature()` now shows the
  effect curve and its curvature with intervals; the "local fragility index"
  and normalised-curvature panels are removed.

## Bug fixes and corrected claims

* The natural indirect and direct effects are now correct when the outcome
  model is not linear in the mediator (for example a logistic outcome).
  Previously the predicted mediator mean was plugged into the outcome model,
  which is exact only for linear outcomes. The outcome prediction is now
  averaged over the mediator's distribution: exactly for a binary mediator and
  by a 15-point quantile rule for a Gaussian mediator.

* `cluster_var` now makes the bootstrap resample whole clusters. Previously it
  only stored cluster sizes and a fixed ICC of 0.05 that was not estimated;
  the ICC field is removed.

* `spline_df` is deprecated and ignored (it never had an effect). Put a spline
  of the treatment in the formulas instead, e.g. `splines::ns(X, df = 3)`.

* `augment()` now returns the same weights used in the balance diagnostics.

* Documentation now states that effects are estimated by g-computation and
  that the inverse probability weights are used for balance diagnostics only.

* Removed unused internal code.

* Repository URLs updated to https://github.com/subirhait/RobustMediate.

# RobustMediate 0.2.0

## Breaking changes

* Removed the mediation ITCV ("medITCV") functions: `sensitivity_meditcv()`,
  `sensitivity_meditcv_profile()`, `plot_meditcv()`, `plot_meditcv_profile()`,
  `fragility_table()`, and their print methods, together with the
  `"meditcv"` and `"meditcv_profile"` options of `plot()`. The quantities they
  reported did not match Frank's (2000) impact threshold for a confounding
  variable, and the "minimum robustness" summary for the indirect effect does
  not hold in general.

* Removed the E-value x rho sensitivity surface, `tipping_table()`, and the
  `evalue_seq` argument of `robustmediate()`. The surface was not computed
  from either the E-value or the sensitivity model of Imai et al. (2010).
  `glance()` now reports `rho_zero_NIE` in place of `evalue_NIE` and
  `rho_min_NIE`.

## New features

* `sensitivity_rho()` implements the sensitivity analysis of Imai, Keele, and
  Yamamoto (2010) for linear mediator and outcome models: the indirect effect
  as a function of the correlation `rho` between the two models' errors, the
  value of `rho` at which it is zero, and the equivalent product of
  variance-explained parameters. `robustmediate()` runs it automatically when
  the models are linear and stores the result in `$sensitivity`.

* `plot_sensitivity()` now plots the indirect effect against `rho`.

## Bug fixes

* `plot_curvature()` failed to draw because the confidence ribbon was applied
  to panels without interval bounds.

* Added the missing `tests/testthat.R` runner, so the test suite now runs
  during `R CMD check`. Fixed a balance test that referred to a non-existent
  field.

# RobustMediate 0.1.1

## Bug fixes

* Fixed bootstrap confidence interval computation. The previous version
  used `stats::update()` to refit models on bootstrap samples, which
  failed silently when bootstrap data frames had duplicate row names
  produced by sampling with replacement. All confidence intervals were
  returning NA. The fix replaces `stats::update()` with direct
  `stats::glm()` calls and resets row names on each bootstrap sample.
  All bootstrap replicates now succeed and confidence intervals are
  correctly computed.

* Fixed balance SMD computation: `max_smd` and `n_above` now handle
  edge cases with near-constant covariates without returning NA.

## Other changes

* Added `n_boot_valid` diagnostic field to `fit$effects` reporting the
  number of successful bootstrap replicates.

* Improved warning messages when bootstrap replicates fail.

* Removed pkgdown URL from DESCRIPTION (was returning 404). Only the
  GitHub repository URL is now listed.

# RobustMediate 0.1.0

* Initial CRAN release.

