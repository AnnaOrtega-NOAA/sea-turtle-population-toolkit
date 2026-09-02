# Sea Turtle Population Toolkit

A manager-facing Shiny workflow for moving from nesting records to population trend estimation and scenario comparison for threats and conservation actions.

## What the app does

1. Upload annual or monthly nesting records.
2. Define population-specific biological parameters.
3. Estimate annual-nester abundance and population trend.
4. Build named threat and conservation scenarios.
5. Compare projected outcomes using matched posterior draws and shared environmental variation.
6. Optionally open **Developer view** for posterior diagnostics and model-detail outputs.

## Methods lineage

The population engine retains the Martin/Siders analytical lineage used by the associated methods work:

- nests are converted to annual nesters using clutch frequency;
- the shared-trend Bayesian state-space model estimates instantaneous trend `U`, process variance `Q`, and observation variance `R`;
- current reproductive-female abundance retains the Siders final-four-year `RI/4` calculation as a reporting quantity;
- fishery Take uses inverse VBGF age, size-dependent maturity, juvenile/adult survival, first-nesting-only ANE, probability female, and remigration interval;
- Take is removed before annual growth;
- projection process error is additive on the abundance scale with SD `sqrt(Q)`;
- conservation Give is an explicit counterfactual extension and is added after the inherited Take/growth step.

The app uses the analytical expected first-nesting probability through time for aggregate-count workflows. This is an expected-value translation of the selected Siders first-nesting-only (`n1`) method, not a claim that every legacy stochastic implementation is reproduced line-for-line.

## Requirements

Install R, JAGS, and the required R packages:

```r
install.packages(c(
  "shiny", "bslib", "dplyr", "tidyr", "ggplot2", "purrr",
  "mvtnorm", "truncnorm", "jagsUI", "MARSS", "cicerone", "car"
))
```

JAGS must also be installed separately on the operating system for the Bayesian model fits.

## Run locally

Open `app.R` in RStudio and click **Run App**, or run:

```r
shiny::runApp("app.R")
```

## Validation

See [`VALIDATION.md`](VALIDATION.md) for the current validation record. The app also contains built-in mathematical self-checks available in Developer view.

Validation checks test implementation consistency and method lineage. They do not establish that any particular population-specific parameter set is empirically adequate.

## Input notes

- Nest counts must be positive where monitoring occurred because the retained trend model is Gaussian on the log-count scale.
- Missing/unmonitored records should be represented as missing rather than silently converted from observed zeroes.
- Threat schedules use affected turtles and a mortality probability; do not pre-multiply counts by mortality.
- Conservation effects are counterfactual benefits relative to what would have happened without the intervention.

## Repository status

This repository is intended as a transparent, reproducible decision-support workflow. The associated manuscript should be cited for full equations, assumptions, scope, and methodological provenance once published.

## License

No open-source license is included yet. Choose and add the appropriate license before inviting unrestricted reuse or redistribution.
# sea-turtle-population-toolkit
