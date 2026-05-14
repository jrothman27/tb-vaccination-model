# TB Vaccination Model

Analysis code for the manuscript:

> Rothman JE, Castro KG, Lopman BA, Gandhi NR, Nelson KN. **Targeted Tuberculosis Vaccination Strategies in the United States: A Modeling Study.** *(Manuscript in preparation.)*
> https://doi.org/10.64898/2026.05.11.26352914

The code in this repository implements a deterministic compartmental transmission model that evaluates targeted tuberculosis vaccination strategies in the United States using an M72/AS01E-like vaccine candidate. The model has four mutually exclusive risk strata (people living with HIV, persons with non-HIV medical comorbidities, non-U.S.-born persons, U.S.-born persons) and compares six vaccination strategies across two parameter scenarios — plausible (VE = 50%, strategy-specific vaccination rates) and optimistic (VE = 70%, vaccination rate = 50%/yr for all strategies). The code was written in R by Jessica Rothman.

## Contents

```
.
├── R/
│   ├── TB_VACCINATION_MODEL_COMPLETE.R       # Base model: ODE system, calibration, shared utilities
│   ├── TB_EXTRACT_CALIBRATION.R              # Calibrated parameter values
│   ├── TB_EXTRACT_STRATEGIES.R               # Per-strategy compartment sizes, summaries, case sources, flows, Mtb infection pool
│   ├── TB_EXTRACT_MASTER.R                   # Merge per-strategy outputs into a single master CSV
│   ├── TB_STRATUM_IMPACT.R                   # Stratum-level cases prevented and NNV
│   ├── TB_TIME_TO_IMPACT.R                   # Annual and cumulative impact over 30 years
│   ├── TB_DIRECT_INDIRECT_PSI_THRESHOLD.R    # Direct/indirect decomposition + psi sensitivity + VE/duration threshold
│   ├── TB_EFFECTIVE_COVERAGE.R               # Effective coverage over time (overall, by stratum, targeted-only)
│   ├── TB_NNV_BENCHMARK_OPTIMISTIC.R         # NNV equilibrium estimates for cross-vaccine benchmarking
│   ├── TB_ONEWAY.R                           # One-way sensitivity sweeps + published NNV comparison
│   ├── TB_LHS.R                              # 500-draw Latin Hypercube Sampling and PRCC
│   ├── TB_CLEAN_ALTERNATIVE_CALIBRATION.R    # Alternative calibration anchored on Aim 1 prevalence
│   ├── TB_DISSERTATION_FIGURES.R             # Generates manuscript and appendix figures
│   └── RUN_ALL.R                             # Master pipeline runner (sources all analyses in order)
├── docs/
│   ├── SCRIPT_OUTPUTS_TABLE.md               # Script-to-output filename map
│   └── CSV_DATA_DICTIONARY.md                # Column-level documentation of every CSV
├── LICENSE                                   # MIT
├── .gitignore
└── README.md
```

## System requirements

- The code is written in R version 4.5.2 using RStudio.
- The following packages are used to run the code: `deSolve`, `lhs`, `sensitivity`, `ggplot2`, `scales`, `patchwork`, `RColorBrewer`, `viridis`, `dplyr`, `tidyr`, `ggtext`.
- The code has been developed and tested on macOS. It should be compatible with Windows, macOS, and Linux operating systems.
- No non-standard hardware is required to run the code. The Latin Hypercube Sampling analysis is the most compute-intensive step and benefits from a multi-core machine, although the code runs serially.

## Installation guide

**Installing the latest version of R**

1. Go to the Comprehensive R Archive Network: <https://cran.r-project.org/>
2. Download the version for your operating system (e.g. Download R for Windows / macOS / Linux).
3. Follow the instructions provided.

This should take about 3 minutes.

**Installing RStudio**

1. Go to the Posit download page: <https://posit.co/downloads/>
2. Download RStudio for your operating system.
3. Follow the instructions provided.

This should take about 2 minutes.

**Installing the R packages**

After installing R and RStudio, install the required packages using the following code:

```
install.packages(c("deSolve", "lhs", "sensitivity", "ggplot2", "scales",
                   "patchwork", "RColorBrewer", "viridis", "dplyr", "tidyr",
                   "ggtext"))
```

This should take about 5 minutes.

## Instructions for running the scripts

Set the working directory to the `R/` folder, then run the full pipeline:

```
setwd("R")
source("RUN_ALL.R")
```

`RUN_ALL.R` sources every analysis script in the correct order and writes all CSV outputs to the working directory. Total runtime is approximately 60 to 240 minutes depending on hardware, dominated by the Latin Hypercube Sampling analysis. To skip individual analyses, comment out the corresponding `source()` calls in `RUN_ALL.R`.

After the pipeline completes, generate the manuscript and appendix figures:

```
source("TB_DISSERTATION_FIGURES.R")
```

To run individual analyses without the full pipeline, each analysis script can be sourced directly. Every script automatically sources `TB_VACCINATION_MODEL_COMPLETE.R` to load the base model and calibration. Most analyses run both scenarios (plausible and optimistic) in a single pass. To restrict to one scenario, comment out the unwanted entry in the `scenarios` list at the top of each script.

## Outputs

The pipeline produces approximately 45 CSV files split into six logical categories:

- **Calibration parameters** — `TB_CALIBRATION_RESULTS.csv`
- **Per-strategy outputs** — compartment sizes, strategy summaries, stratum-level impact, case sources, *Mtb* infection pool, annual flows, master merge (7 CSVs per scenario)
- **Sensitivity analyses** — one-way parameter sweeps, tornado summaries, LHS uncertainty intervals, PRCC, published NNV comparison
- **Threshold analyses** — full VE × duration grid, minimum-VE table
- **Coverage analyses** — psi dose-response, effective coverage (overall, by stratum, targeted-only)
- **Time-varying outputs** — time-to-impact trajectories, direct/indirect decomposition, NNV benchmark

Column-level documentation for every CSV is in `docs/CSV_DATA_DICTIONARY.md`. A script-to-output filename map is in `docs/SCRIPT_OUTPUTS_TABLE.md`.

The figures script generates all manuscript figures (Figures 1–3) and appendix figures (Figures S1–S5) as both `.pdf` and `.png` at 300 dpi.

## Data

TB case counts are from the U.S. National Tuberculosis Surveillance System (NTSS), 2024 (10,388 reported cases). The model is calibrated by fitting the transmission rate (β) and four stratum-specific *Mtb* infection prevalences among new entrants to match 2024 NTSS stratum-specific incidence targets. Stratum-specific population sizes are from the CDC HIV Surveillance Report (PLWH), CDC National Diabetes Statistics Report and related sources (medical comorbidities), American Community Survey (non-U.S.-born), and U.S. Census (U.S.-born). All parameter values and sources are documented in Tables S1–S3 of the manuscript appendix. All inputs needed to reproduce the results are hard-coded in `TB_VACCINATION_MODEL_COMPLETE.R`.

## License

Code is released under the MIT License (see `LICENSE`).

## Contact

Jessica E. Rothman — jessica.e.rothman@gmail.com
