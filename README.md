# 🐢 Sea Turtle Population Toolkit

<p align="left">
  <img src="https://img.shields.io/badge/R-%3E%3D4.0-276DC3.svg?logo=r" alt="R version" />
  <img src="https://img.shields.io/badge/Shiny-Interactive%20App-0072B2.svg?logo=rstudio" alt="Shiny" />
  <img src="https://img.shields.io/badge/Backend-CmdStanR-CC0000.svg" alt="CmdStanR" />
  <img src="https://img.shields.io/badge/Status-Early%20Public%20Release-2E7D32.svg" alt="Status" />
</p>

> **A decision-support application for translating nesting-beach monitoring records into unified population trajectories and evaluating management scenarios for threats and conservation actions.**

Designed for conservation practitioners, managers, and researchers who need an evidence-based path from incomplete field monitoring data to population-level viability without requiring custom Bayesian MCMC coding.

---

## 🧭 The 5-Stage Workflow

1. **Add Nesting Data**  
   Upload annual or monthly nesting records across single or multiple index beaches. The built-in Bayesian Fourier model reconstructs unmonitored seasonal gaps while preserving observed counts.

2. **Describe Population Biology**  
   Specify population-specific demographic parameters (clutch frequency, remigration interval, somatic growth, and stage survival) rather than defaulting to generic global averages.

3. **Estimate Status & Trend**  
   Fit a unified Bayesian state-space model via **CmdStanR** to estimate long-term exponential growth rate (*r*), annual process error (*σ*<sub>proc</sub>), and site-specific scaling factors.

4. **Build Scenarios (Take & Give)**  
   Define threat schedules (mortality across size classes) and conservation actions (nest protection, headstarting, adult mortality prevention) standardized into a shared metric: **Annual Nester Equivalents (ANE)**.

5. **Project & Compare Outcomes**  
   Run paired stochastic 30-year projections. Alternative management actions are evaluated against the same environmental noise and posterior draws to isolate real policy differences from simulation noise.

---

## ⚖️ Demographic Standardization: Annual Nester Equivalents (ANE)

Monitoring counts breeding females, but threats and conservation actions often impact juveniles, subadults, or eggs. Directly comparing raw counts across life stages obscures survival to maturity, somatic growth delays, and remigration intervals. The toolkit standardizes all impacts into **Annual Nester Equivalents (ANE)**:

* **Threats (Take):** An interaction at length *L* is projected through somatic growth to determine the timing and probability of first nesting, discounted by cumulative survival, sex ratio, remigration interval, and interaction mortality.
* **Conservation (Give):** Credited strictly as the *incremental benefit* relative to a no-intervention counterfactual. Additional hatchlings, headstarted yearlings, or saved adults follow the same demographic engine before contributing to future nesting abundance.

<details>
<summary><b>🔬 View Mathematical Formulations & Projection Equations</b> (Click to expand)</summary>

<br>

### 1. Somatic Growth (von Bertalanffy)
$$L(A) = L_\infty \left[ 1 - \exp\{-k(A - t_0)\} \right]$$

### 2. Transition Toward First Nesting
$$p(T_{i,t}) = \begin{cases} 1, & L_{i,t} \ge 0.99 L_\infty \\ \left[ 1 + \exp\left( -\frac{L_{i,t} - L_{\text{mat}}}{\sigma_{\text{mat}}} \right) \right]^{-1}, & L_{i,t} < 0.99 L_\infty \end{cases}$$

### 3. First-Nesting Distribution (Discrete Hazard)
Calculates the marginal probability that first nesting occurs in year $t$ without Monte Carlo sampling error:
$$p(FN_{i,t}) = p(T_{i,t}) \prod_{y=1}^{t-1} \left[ 1 - p(T_{i,y}) \right]$$

### 4. Cumulative Survival
$$\Phi_{i,t} = \left[ 1 - p(T_{i,t}) \right] \Phi_{\text{juv}} + p(T_{i,t}) \Phi_{\text{adult}}, \qquad S_{i,t} = \prod_{y=1}^t \Phi_{i,y}$$

### 5. Annual Nester Equivalent (ANE) Standardization
* **Threat (Take):**
  $$ANE_{i,t}^T = S_{i,t} \cdot p(FN_{i,t}) \cdot p(F) \cdot \left(\frac{1}{RI}\right) \cdot D_i$$
* **Conservation (Give):**
  $$ANE_{a,t}^G = \Delta n_a \cdot S_{a,t} \cdot p(FN_{a,t}) \cdot p(F) \cdot \left(\frac{1}{RI}\right)$$

### 6. Recursive Population Projection
$$N_{t+1} = \max\left[0, \; \max\left(0, N_t - ANE_t^T\right) \exp\left(r + \sigma_{\text{proc}} z_t\right) + ANE_{t+1}^G\right]$$

where $z_t \sim \text{Normal}(0, 1)$ represents annual environmental process variation. Take is deducted prior to annual growth, while matured Give is added after the transition.

</details>

---

## 📋 Data File Formats

Upload data using clean CSV files formatted as follows:

### 1. Nesting Monitoring Data
```csv
Year,Month,Site,Count,Monitored
2015,5,Beach_A,42,TRUE
2015,6,Beach_A,,FALSE
2015,7,Beach_A,0,TRUE
```

> [!NOTE]
> Distinguish periods with no survey effort (`Count = blank, Monitored = FALSE`) from confirmed zero nesting (`Count = 0, Monitored = TRUE`).

### 2. Threat Schedules (Take)
```csv
year,turtles,median_cm,mortality
2028,25,85,0.35
2029,15,142,1.00
```

### 3. Conservation Schedules (Give)
```csv
Year,Amount
2028,500
2029,500
```
*(Where `Amount` represents clutches protected, yearlings released, or adult mortalities averted).*

---

## ⚙️ Installation & Requirements

The toolkit runs in **R** ($\ge 4.0$) and uses **CmdStanR** for Bayesian model fitting.

### 1. Install C++ Toolchain
Stan compiles models to native C++. Ensure a C++ compiler is available:
* **Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version.
* **macOS:** Run `xcode-select --install` in the terminal.
* **Linux:** Install `build-essential` and `g++`.

### 2. Install CmdStanR & CmdStan Backend
```r
# Install cmdstanr from the Stan R-universe
install.packages("cmdstanr", repos = c("[https://stan-dev.r-universe.dev](https://stan-dev.r-universe.dev)", getOption("repos")))

# Verify toolchain and install CmdStan
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
```

### 3. Install R Packages
```r
install.packages(c(
  "shiny", "bslib", "dplyr", "tidyr", 
  "ggplot2", "purrr", "mvtnorm", "truncnorm", 
  "posterior", "coda", "MARSS", "car"
))
```

---

## 🏃 Launching the Application

Clone the repository and run the app from your local terminal or R console:

```bash
git clone [https://github.com/AnnaOrtega-NOAA/sea-turtle-population-toolkit.git](https://github.com/AnnaOrtega-NOAA/sea-turtle-population-toolkit.git)
cd sea-turtle-population-toolkit
```

In R or RStudio:
```r
shiny::runApp()
```

> [!TIP]
> **First-Time Walkthrough:**  
> 1. Launch the app and click **Load Synthetic Dataset**.  
> 2. Step through the guided tabs: estimate trend, build an illustrative threat scenario and a conservation scenario, and run the comparative projection.  
> 3. Toggle **Developer View** in the top navigation to inspect MCMC trace plots, posterior chains, and convergence diagnostics.

---

## 🔬 Computational & Mathematical Verification

The analytical engine has undergone extensive numerical validation:

* **Analytical Consistency:** Confirmed machine-precision numerical agreement against historical assessment implementations.
* **Biological Symmetry:** Enforced structural conservation laws (e.g., preventing one adult mortality exactly offsets one adult lethal take; cohort nesting probabilities sum to 1.0).
* **Monitoring Stress Tests:** Demonstrated that Bayesian trend estimation remains robust under up to 60% random monthly survey loss, while highlighting that multi-year monitoring gaps primarily degrade terminal abundance estimates.

---

## 📚 References

The analytical framework synthesizes NOAA Pacific Islands Fisheries Science Center assessments and recent methodological advancements:

* **Martin, S. L., Siders, Z., Eguchi, T., Langseth, B. J., Yau, A., Baker, J. D., Ahrens, R. N. M., & Jones, T. T. (2020a).** *Assessing the Population-level Impacts of North Pacific Loggerhead and Western Pacific Leatherback Turtle Interactions in the Hawaii-based Shallow-set Longline Fishery*. NOAA Technical Memorandum NMFS-PIFSC-95, 137 p. [doi:10.25923/ydp1-f891](https://doi.org/10.25923/ydp1-f891)
* **Martin, S. L., Siders, Z., Eguchi, T., Langseth, B., Yau, A., Baker, J., Ahrens, R., & Jones, T. T. (2020b).** *Update to Assessing the Population-level Impacts of North Pacific Loggerhead and Western Pacific Leatherback Turtle Interactions: Inclusion of the Hawaii-based Deep-set and American Samoa-based Longline Fisheries*. NOAA Technical Memorandum NMFS-PIFSC-101, 67 p. [doi:10.25923/pnf2-2q77](https://doi.org/10.25923/pnf2-2q77)
* **Siders, Z. A., Martin, S. L., Ahrens, R. N. M., & Jones, T. T. (2023).** *Update incorporating uncertainty in maturation and fishery takes into population-level impacts of western Pacific leatherback sea turtles*. PIFSC Internal Report IR-23-03, 29 p.
* **Siders, Z. A., Martin, S. L., Ahrens, R. N. M., & Jones, T. T. (2025).** *Update to NOAA Technical Memorandum NMFS-PIFSC-101: Incorporating Uncertainty in Maturation and Recent Fishery Takes into the Population-Level Impacts of Western Pacific Leatherback Sea Turtles*. NOAA Technical Memorandum TM-PIFSC-184, 43 p. [doi:10.25923/72944](https://doi.org/10.25923/72944)
* **Ortega, A. A., Siders, Z. A., Eguchi, T., Ahrens, R. N. M., Jones, T. T., & Martin, S. L. (In review).** *A common demographic currency for linking anthropogenic impacts to population status in sea turtles*. Conservation Biology.

---

---

## ⚖️ License

Software authored by U.S. Government employees is in the public domain within the United States pursuant to 17 U.S.C. § 105. 

For international use and contributions authored by non-federal partners (e.g., University of Hawaiʻi / CIMAR), this project is distributed under the **MIT License** (or dedicated to the public domain via **CC0 1.0 Universal**). See [LICENSE.md](LICENSE.md) for full terms and the Department of Commerce copyright addendum.

---

## 🏛️ NOAA & Department of Commerce Disclaimer

This repository is a scientific product and is not an official communication of the National Oceanic and Atmospheric Administration (NOAA), the National Marine Fisheries Service (NMFS), or the United States Department of Commerce (DOC). 

All NOAA GitHub project code is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

The scientific results, models, and conclusions presented in this software are intended to support decision analysis and do not represent formal regulatory determinations, Endangered Species Act (ESA) biological opinions, or official agency policy.
