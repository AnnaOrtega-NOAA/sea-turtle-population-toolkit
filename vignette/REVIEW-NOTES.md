# Sea Turtle Population Toolkit — redesign review

Scientific sources inspected:
- Methods_ConsBio_20260831.pdf
- MethodsScript_ConsBio_20260831(8).R
- app(20260902-172448).R

The redesign changes the Quarto teaching site, not the Shiny population engine.

## Verification completed
- Quarto 1.8.24: `quarto preview` started and watched the project.
- Quarto 1.8.24: `quarto render` completed for all three pages.
- 19 automated scientific and DOM-interaction test groups pass.
- Five paper-reported cumulative ANE values reproduced at their printed precision.
- Every data-entry branch, concept link, disclosure, size preset, conservation input group, and matched-future control exercised by automated DOM tests.
- Initial browser review: convergence page, data branches, and adult size explorer inspected.
- Initial browser review detected Quarto's `.columns` collision. Renamed the layout class to `.science-columns` and rerendered.
- Internal file links and MathML equation output checked after rendering.

## Remaining review gate
The browser's security policy blocked further preview access after the layout fix. Final desktop/mobile visual QA and a complete live-browser interaction pass were therefore not completed. The site has not been deployed. Open `quarto preview` locally to review the final layout before publication.

The supplied app's 25/25 mathematical validation record is attributed to its own header. R/JAGS tests were not rerun for this website redesign. The site's 19 automated checks are separate from that record.

## Use the files
Overlay `_quarto.yml`, `index.qmd`, `walkthrough.qmd`, `styles.css`, `assets/`, and `guide/developer.qmd` in the existing Quarto website. Keep the existing Shiny app and other repository files.

Run `quarto preview` for a local review or `quarto render` to rebuild `_site/`. The site is static; no R or JAGS fitting is performed. `npm install` followed by `npm test` runs the teaching tests. Node dependencies are only needed for tests and the optional rendered-output preview, not Quarto rendering.

The package includes a rendered `_site/` copy. Serve it through a local HTTP server; opening module-based interactions directly with file:// is not supported by browsers.
