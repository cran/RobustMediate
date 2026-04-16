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

