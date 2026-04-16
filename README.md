# RobustMediate

**Robust causal mediation analysis** with embedded diagnostics, 
dose-response curves, pathway-specific sensitivity (medITCV), 
and a novel bivariate sensitivity contour.

## What it does

| Function | What it gives you |
|---|---|
| `robustmediate()` | Fit treatment / mediator / outcome models, compute IPW weights, NDE/NIE/TE curves with bootstrap CIs, and the full sensitivity surface in one call |
| `plot_balance()` | Dual love plot: covariate balance before/after weighting for **both** pathways simultaneously |
| `plot_mediation()` | Dose-response curves of NDE, NIE, TE with pointwise confidence bands |
| `plot_sensitivity()` | Novel 2-D robustness map: E-value x Imai rho — does not exist elsewhere in R |
| `sensitivity_meditcv()` | Pathway-specific mediation ITCV (medITCV) for a-path and b-path |
| `plot_meditcv()` | Robustness corridor plot for each pathway |
| `sensitivity_meditcv_profile()` | Minimum robustness principle + bottleneck identification |
| `plot_meditcv_profile()` | Fragility profile as confounding impact increases |
| `fragility_table()` | Publication-ready pathway decomposition table |
| `diagnose()` | Formatted report with a paste-ready Results paragraph |

## Installation
```r
# Development version from GitHub
# install.packages("pak")
pak::pkg_install("causalfragility-lab/RobustMediate")
```

## Quick start
```r
library(RobustMediate)

fit <- robustmediate(
  treatment_formula = X ~ Z1 + Z2,
  mediator_formula  = M ~ X + Z1 + Z2,
  outcome_formula   = Y ~ X + M + Z1 + Z2,
  data = mydata,
  R    = 500
)

plot_balance(fit)                        # love plot
plot_mediation(fit)                      # NDE / NIE dose-response curve
plot_sensitivity(fit)                    # E-value x rho contour
plot(fit, type = "meditcv")             # medITCV robustness corridor
plot(fit, type = "meditcv_profile")     # fragility profile
fragility_table(fit)                     # pathway decomposition
diagnose(fit)                            # paste into Results section
```

## Why this package?

- `EValue` plots E-values only
- `mediation` plots rho sensitivity only  
- `cobalt` / `WeightIt` do love plots for treatment only

**RobustMediate** combines all three into one coherent workflow 
tailored to continuous-treatment mediation, and adds:

- The joint E-value x rho contour that exists nowhere else in R
- Pathway-specific medITCV (mediation ITCV) extending Frank (2000) to mediation
- Minimum robustness principle and bottleneck identification for indirect effects

## References

- Frank, K. A. (2000). Impact of a confounding variable on a regression coefficient. *Sociological Methods & Research*, 29(2), 147-194.
- VanderWeele, T. J. & Ding, P. (2017). Sensitivity analysis in observational research: Introducing the E-value. *Annals of Internal Medicine*, 167(4), 268-274.
- Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference and sensitivity analysis for causal mediation effects. *Statistical Science*, 25(1), 51-71.

## Contributing

Bug reports and feature requests via [GitHub Issues](https://github.com/causalfragility-lab/RobustMediate/issues).

## License

MIT

