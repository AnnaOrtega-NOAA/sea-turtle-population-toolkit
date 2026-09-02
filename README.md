# Sea Turtle Population Toolkit

A Shiny decision-support toolkit for turning sea turtle nesting records into population trends and scenario comparisons for threats and conservation actions.

The toolkit is designed for conservation practitioners who need a transparent path from monitoring data to population-level inference without requiring advanced modeling expertise.

---

## Overview

Sea turtle monitoring programs often collect nesting data over many years, but those records can contain missing months, multiple beaches, changing survey effort, and substantial uncertainty.

At the same time, management questions are rarely limited to whether a population is increasing or decreasing. Practitioners may also need to ask:

- What could happen if a threat continues?
- How much difference could reducing mortality make?
- When would benefits from nest protection or headstarting become visible?
- How do several management scenarios compare?
- How much uncertainty comes from the monitoring data or demographic assumptions?

The Sea Turtle Population Toolkit brings these steps into one workflow.

Users can move from nesting records to population estimation, construct threat or conservation scenarios, and compare projected outcomes using a shared demographic currency.

---

## What the toolkit does

The application guides users through five main stages.

### 1. Add nesting data

Upload annual or monthly nesting records.

The toolkit can:

- use annual nest counts directly;
- accommodate multiple nesting sites;
- identify periods without monitoring;
- reconstruct missing monthly observations using seasonal Fourier imputation;
- distinguish missing monitoring from genuine observed counts.

### 2. Describe the population

Enter population-specific biological parameters such as:

- clutch frequency;
- remigration interval;
- proportion female;
- von Bertalanffy growth parameters;
- size at maturity;
- maturity-ogive width;
- juvenile survival;
- adult survival.

The application does not automatically substitute global species averages for population-specific parameters.

### 3. Estimate population status

The toolkit converts nests to annual nesting females and estimates:

- historical annual-nester abundance;
- long-term instantaneous population trend;
- process variance;
- observation variance;
- regional population trajectories;
- current reproductive-female abundance summaries.

An optional developer view exposes additional diagnostics and posterior outputs.

### 4. Build scenarios

Users can save multiple named scenarios representing:

- threats;
- conservation actions;
- combined portfolios.

Threat examples include fisheries interactions or other sources of mortality.

Conservation examples can include:

- nest protection;
- headstarting;
- prevention of adult mortality.

Scenarios are retained independently so alternative assumptions can be saved and compared.

### 5. Compare projected outcomes

Saved scenarios are projected against the same underlying population state.

Scenario comparisons share:

- initial posterior draws;
- annual trend and process-variance draws;
- environmental stochasticity.

This matched design reduces Monte Carlo noise when comparing management alternatives.

Outputs include:

- projected annual nesting-female abundance;
- uncertainty intervals;
- scenario summaries;
- downloadable simulation draws;
- run settings;
- mathematical diagnostics.

---

## Scientific framework

The toolkit retains the mathematical lineage of population assessments developed by Martin and subsequently extended by Siders for sea turtle population-impact analyses.

The retained components include:

- conversion of nesting records to annual nesting females;
- Bayesian state-space estimation of population trend;
- process and observation variance;
- annual-nester population projections;
- remigration-interval abundance summaries;
- translation of affected turtles into annual nester equivalents;
- first-nesting-only demographic accounting;
- fishery Take applied before population growth;
- paired dynamic posterior draws of population trend and process variance.

The application also introduces explicitly labeled extensions for conservation actions, scenario comparison, input uncertainty, and decision-support presentation.

---

## Annual nester equivalents

Threats affect turtles at different ages and sizes, while the monitored population is usually expressed in terms of nesting females.

The toolkit therefore translates affected turtles into **annual nester equivalents (ANE)**.

For a turtle entering the calculation at a given size, the expected contribution to future annual nesting abundance depends on:

- growth and age;
- probability of maturity;
- probability of first nesting;
- survival;
- probability of being female;
- remigration interval;
- mortality associated with the threat.

The probability of first nesting in year \(j\) is represented as:

```math
p(FN_j)
=
p(M_j)
\prod_{h<j}
\left[1-p(M_h)\right]
```

This assigns each turtle to its expected first nesting event only once.

The resulting threat contribution is standardized to annual nesting females.

---

## Population projection

The retained Martin/Siders projection structure is:

```math
N_{t+1}
=
\max\left[
0,\;
\max(0,N_t-T_t)e^{U_t}
+
\sqrt{Q_t}z_t
+
G_t
\right]
```

where:

- \(N_t\) is annual nesting-female abundance;
- \(T_t\) is threat-related ANE loss;
- \(U_t\) is instantaneous population trend;
- \(Q_t\) is process variance;
- \(z_t\) is a standard-Normal process innovation;
- \(G_t\) is conservation-related ANE gain.

Threat-related Take is removed **before population growth**.

Conservation Give is an explicit extension and is added **after the retained Take-growth-process transition**.

Negative projected abundance is bounded at zero as a numerical safeguard.

---

## Conservation Give

Positive conservation effects are represented as **counterfactual gains**.

The toolkit does not automatically credit every turtle, egg, hatchling, or adult associated with a conservation program.

Instead, the relevant quantity is the additional contribution attributable to the intervention relative to what would have happened without it.

### Nest protection

Benefit is based on the additional hatchlings produced relative to the unprotected counterfactual.

### Headstarting

Benefit is based on the **net additional turtles alive because of the intervention**, rather than simply the number released.

### Adult protection

Benefit is based on adult deaths prevented relative to the counterfactual.

These additional cohorts are propagated forward through survival, maturation, first nesting, sex ratio, and remigration rather than receiving an immediate adult-equivalent credit.

---

## Monitoring-data gaps

Monthly nesting records can contain periods when surveys did not occur.

The toolkit distinguishes:

- a genuine observed count;
- an observed zero;
- a period with no monitoring.

Missing monitoring can be reconstructed using the seasonal imputation model.

A genuine observed zero is not silently treated as missing because the retained Gaussian log-count model cannot accommodate zero directly.

Users should therefore preserve the distinction between:

```text
Count = blank, Monitored = FALSE
```

and:

```text
Count = 0, Monitored = TRUE
```

The latter requires a different observation-model formulation.

---

## Expected input format

### Nesting data

CSV columns:

```text
Year,Month,Site,Count,Monitored
```

For annual records, `Month` can be left blank.

`Count` refers to nests, not turtles.

### Threat schedules

CSV columns:

```text
year,turtles,median_cm,mortality
```

where:

- `turtles` = number affected during that year;
- `median_cm` = representative carapace length;
- `mortality` = probability of death per affected turtle.

For example:

```text
2030,25,85,0.35
```

means 25 turtles affected, with representative length 85 cm and mortality probability 0.35.

Do not pre-multiply turtle numbers by mortality.

### Conservation schedules

Conservation uploads use:

```text
Year,Amount
```

The meaning of `Amount` depends on the selected action.

Examples:

- nest protection: nests receiving protection;
- headstarting: net additional turtles alive because of the intervention;
- adult protection: adult deaths prevented.

The interface collects the biological assumptions needed to translate those quantities into delayed ANE contributions.

---

## Guided and developer views

The default interface is designed for practitioners.

It emphasizes:

- plain-language questions;
- progressively revealed inputs;
- compact summaries;
- scenario-based decision making;
- minimal statistical jargon.

A **Developer view** can be enabled for technical inspection.

Developer outputs include additional information on:

- posterior distributions;
- convergence diagnostics;
- process and observation variance;
- model comparison;
- retained posterior chains;
- mathematical checks.

The scientific engine is the same in both views. Developer view changes presentation, not the underlying calculations.

---

## Mathematical validation

The current methods-aligned implementation includes regression and identity checks covering:

- projection-equation identity;
- zero-Take equivalence;
- Take-before-growth timing;
- zero interactions;
- zero mortality;
- stochastic reproducibility;
- first-nesting accounting;
- ANE response to interaction number and mortality;
- zero conservation action;
- linear conservation-action scaling;
- conservation timing;
- emergence-success denominator equivalence;
- adult Give/Take counterfactual symmetry;
- Siders RI/4 abundance identity;
- matched posterior \(U/Q\) rows.

During development, the mathematical engine passed:

```text
565 external reference comparisons
25 built-in mathematical checks
```

with a maximum reported reference difference of approximately:

```math
1.39 \times 10^{-17}
```

The guided workflow was additionally exercised using annual and monthly examples, JAGS/MARSS fitting, scenario creation, comparison, downloads, and mathematical diagnostics.

These checks evaluate software implementation and mathematical consistency.

They do **not** demonstrate that a particular set of biological parameters is appropriate for a particular population.

---

## Important interpretation note

This software should be described as **method-lineage aligned**, rather than as a literal reproduction of every historical Martin or Siders implementation.

In particular, aggregate fishery ANE in the interactive framework uses the analytical expected value of the selected Siders first-nesting-only formulation.

The application retains the corresponding:

- growth model;
- maturation structure;
- survival calculation;
- female probability;
- remigration standardization;
- Take timing;
- dynamic paired \(U/Q\) projection structure.

Conservation Give, user-facing scenario comparison, and selected uncertainty options are explicit extensions.

---

## Installation

The application requires R and JAGS.

### R

A recent version of R is recommended.

### JAGS

Install JAGS separately before running the application.

JAGS is available from:

https://mcmc-jags.sourceforge.io/

### R packages

Install the required packages in R:

```r
install.packages(c(
  "shiny",
  "bslib",
  "dplyr",
  "tidyr",
  "ggplot2",
  "purrr",
  "mvtnorm",
  "truncnorm",
  "jagsUI",
  "MARSS",
  "car"
))
```

---

## Run the application

Clone the repository:

```bash
git clone https://github.com/AnnaOrtega-NOAA/sea-turtle-population-toolkit.git
```

Move into the repository:

```bash
cd sea-turtle-population-toolkit
```

Then launch R and run:

```r
shiny::runApp("app.R")
```

Alternatively, open `app.R` in RStudio and click **Run App**.

---

## Recommended first run

For a first run:

1. Launch the app.
2. Choose the synthetic example dataset.
3. Work through the guided interface.
4. Estimate the population.
5. Create one threat scenario.
6. Create one conservation scenario.
7. Save both scenarios.
8. Compare them with status quo.
9. Enable Developer view if you want to inspect diagnostics.

Once the example workflow is familiar, repeat the process using your own nesting records and population-specific demographic parameters.

---

## Reproducibility

For analytical or management use, retain:

- the original monitoring CSV;
- biological parameter values and their sources;
- scenario schedules;
- uncertainty assumptions;
- downloaded run settings;
- the software version or Git commit used for the analysis.

For formal analyses, sufficiently large simulation and MCMC settings should be used to obtain stable estimates.

Fast settings intended for interactive exploration should not automatically be treated as assessment-quality settings.

---

## Scope

The toolkit is intended to support transparent exploration of population-level consequences of threats and conservation actions.

It is **not** intended to:

- provide universal demographic parameters for sea turtle species;
- replace population-specific biological review;
- establish species-specific extinction thresholds;
- determine regulatory significance automatically;
- substitute for expert assessment or management review;
- infer intervention cost-effectiveness.

Outputs should be interpreted in the context of the monitoring data, demographic evidence, and assumptions supplied by the user.

---

## Method lineage

The analytical framework builds on work including:

- Martin, S. L. et al. (2020). Population-level assessment methods for North Pacific loggerhead and western Pacific leatherback turtle interactions. NOAA Pacific Islands Fisheries Science Center.
- Siders, Z. A. et al. (2023). Update incorporating uncertainty in maturation and fishery takes. PIFSC Internal Report IR-23-03.
- Siders, Z. A., Martin, S. L., Ahrens, R. N. M., & Jones, T. T. (2025). Update incorporating uncertainty in maturation and recent fishery takes into population-level impacts of western Pacific leatherback sea turtles. NOAA Technical Memorandum TM-PIFSC-184.

Please consult the associated reports and manuscripts for methodological detail and appropriate interpretation.

---

## Citation

If you use this toolkit in research, assessment, or management work, please cite the software version used.

Citation information is provided in:

```text
CITATION.cff
```

A manuscript describing the generalized common-currency framework and conservation extensions is in development.

---

## Development status

**Current release stage: early public release / technical review**

The core mathematical engine has undergone extensive implementation testing, but the software remains under active development.

Feedback is particularly useful on:

- usability for conservation practitioners;
- population-specific parameter requirements;
- monitoring-data edge cases;
- diagnostic presentation;
- threat and conservation scenario definitions;
- reproducibility and reporting outputs.

---

## Contributing

Issues and suggestions are welcome.

When reporting a problem, please include:

- the software version or Git commit;
- the stage of the workflow where the problem occurred;
- the relevant error message;
- a minimal reproducible example when possible.

Do not upload sensitive or restricted monitoring data to a public GitHub issue.

---

## License

A software license has not yet been assigned.

Please check the repository before reuse or redistribution, as licensing information may change before the first formal release.

---

## Disclaimer

This software is a scientific decision-support and research tool.

Outputs depend on the quality and appropriateness of the monitoring data, biological parameters, model assumptions, and scenarios supplied by the user.

The toolkit does not independently validate population-specific demographic assumptions or determine management decisions.
