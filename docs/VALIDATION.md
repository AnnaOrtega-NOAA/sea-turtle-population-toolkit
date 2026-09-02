# Guided UI revision — 2 September 2026

Use `app_guided_ui.R` as the app entry point. Open it in RStudio and click **Run App**. It is a self-contained R file; no companion data folder is required. The CSVs in `turtle_pva_upload_templates.zip` remain compatible.

## What changed

- Five steps: add nesting data, set population biology, estimate the population, build scenarios, compare results.
- Persistent form controls and guarded navigation. Switching pages does not recreate the scenario forms.
- Populated examples beside each upload. Annual and monthly monitoring examples can also be loaded directly. Threat and action examples use the selected forecast years.
- Immediate upload validation, a schedule preview, and a plain-language review of amounts, size, mortality, sex assumptions and conservation counterfactuals.
- Named scenarios can be saved, edited, copied and removed. Separate conservation uploads retain their own files. Saved scenarios retain CSV contents rather than relying on temporary upload paths.
- One shared forecast configuration and one comparison run. Scenarios share initial posterior draws, paired annual trend/variance draws and environmental shocks.
- Input uncertainty is beside the relevant scenario inputs; methodological options have a local disclosure. There is no global advanced checkbox.
- Results are flagged when saved settings change. Unsaved edits to a selected scenario are identified before comparison. Stale result downloads are blocked.
- A comparison table and downloadable plot, summary, simulations and human-readable run settings.
- Sensitivity exploration uses named scenario copies in the same comparison engine instead of the former separate sensitivity plot workflow.

## Mathematical scope

The growth, first-nesting, survival, threat, conservation and population-projection kernels retain the methods-aligned equations. Scenario effect construction was extracted from the existing main engine and receives saved input snapshots.

Integration testing found and repaired one additional pre-existing error: the multi-site monthly Fourier JAGS model string had a redundant closing brace. The repair changes the model string's syntax, not its equations. Directly loaded monthly examples now enter the same imputation pathway as uploaded monthly records.

The interface requires explicit rows for every forecast year in both threat and action CSVs. Action uploads with fractional years, negative amounts, missing years or out-of-range years are rejected before the existing alignment helper is called.

## Validation completed

- R parsing and UI construction; no duplicate HTML element IDs.
- 565 reference comparisons/checks and all 25 built-in mathematical checks passed. Maximum difference in the external reference comparisons: 1.39e-17.
- Reactive server tests exercised named scenario creation, independent saved records, uploaded-content retention after removal of the temporary file, shared comparison draws, threat/action direction, unsaved edits, invalid schedule years and stale-result detection.
- Real annual JAGS/MARSS fitting on the 2006–2025 synthetic monitoring example; subsequent projections included a threat and all three conservation action uploads.
- Real monthly Fourier imputation and JAGS/MARSS fitting on the 480-row monthly example, with April season start and the median imputation branch.
- A desktop Chromium walkthrough exercised example loading, biology selection, real historical fitting, threat and nest-protection CSV uploads, named scenario saving, comparison, summary download and mathematical diagnostics. No visible Shiny output errors remained at completion.
- Screenshots of the principal workflow screens were inspected for layout and readability.

`validate_guided_ui.R` provides the repeatable reactive/format tests. Run it with:

```r
# From a terminal in the directory containing the downloaded files:
Rscript validate_guided_ui.R app_guided_ui.R
```

It requires the app dependencies plus `testthat`. It uses a controlled posterior fixture for UI orchestration and does not refit JAGS. The separate `validate_methods_alignment.R` file previously supplied remains the manuscript-reference check suite.

Validation used R 4.3.3 and JAGS 4.3.2. Demonstration fits and small projection runs test functionality; their convergence and Monte Carlo precision are not assessment recommendations. The long legacy MCMC configurations, all uncertainty combinations, and every browser/device were not exercised in this UI revision.

## Input and storage boundaries

Threat CSVs use `year,turtles,median_cm,mortality`. Conservation CSVs use `Year,Amount`, with action-specific meaning shown beside the upload. Sex and other biological assumptions remain shared controls rather than new CSV columns. No unapproved cohort-schema redesign is included.

Scenarios persist within the running app session. The run-settings text download records assumptions and schedules; it is not a session-import file. Closing/restarting the app clears session state.

All embedded examples are synthetic. Model validation does not establish the suitability of biological assumptions for a real population.
