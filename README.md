# RobustMediate

**Causal mediation analysis** for continuous treatments with embedded
balance diagnostics, dose-response curves, and sensitivity analysis for
mediator-outcome confounding.

## What it does

| Function | What it gives you |
|---|---|
| `robustmediate()` | Fit treatment / mediator / outcome models, NDE/NIE/TE curves by g-computation with bootstrap CIs, and the sensitivity analysis in one call |
| `plot_balance()` | Dual love plot: covariate balance before/after weighting for **both** pathways simultaneously |
| `plot_mediation()` | Dose-response curves of NDE, NIE, TE with pointwise confidence bands |
| `sensitivity_paths()` | Robustness values for each path and for the indirect effect under joint-significance and Sobel tests (Cinelli & Hazlett, 2020) |
| `sensitivity_rho()` | Indirect effect as a function of the mediator-outcome error correlation `rho` (Imai, Keele & Yamamoto, 2010) |
| `plot_sensitivity()` | Plot of the indirect effect against `rho` |
| `diagnose()` | Formatted report with a paste-ready Results paragraph |

## Installation
```r
# Development version from GitHub
# install.packages("pak")
pak::pkg_install("subirhait/RobustMediate")
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
plot_sensitivity(fit)                    # indirect effect vs rho
diagnose(fit)                            # paste into Results section
```

## Why this package?

`cobalt` and `WeightIt` give balance diagnostics for the treatment model only,
and `mediation` targets binary or categorical treatments. **RobustMediate**
puts balance checks for both pathways, dose-response curves for a continuous
treatment, and a sensitivity analysis into one workflow.

## References

- Imai, K., Keele, L., & Yamamoto, T. (2010). Identification, inference and sensitivity analysis for causal mediation effects. *Statistical Science*, 25(1), 51-71.

## Contributing

Bug reports and feature requests via [GitHub Issues](https://github.com/subirhait/RobustMediate/issues).

## License

MIT

