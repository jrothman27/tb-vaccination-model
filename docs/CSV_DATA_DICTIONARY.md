# CSV Output Data Dictionary

This document describes the contents of every CSV file produced by the analysis pipeline. Each section lists the file's columns, what each column represents, units where applicable, and the typical row count.

For each CSV, "Strategy" or "Scenario" refers to the vaccination strategy being evaluated (or "Baseline" for the no-vaccination reference). "Stratum" refers to one of the four risk strata: HIV (PLWH), Medical (medical comorbidities), FB or NUSB (non-U.S.-born), USB (U.S.-born).

---

## Calibration outputs

### `TB_CALIBRATION_RESULTS.csv`

Calibrated parameter values from the base model. Scenario-independent.

**Source:** `TB_EXTRACT_CALIBRATION.R` (also produced as a byproduct of `extract_all_model_outputs.R`)

| Column | Description |
|---|---|
| `Parameter` | Parameter name (e.g., `beta`, `mtb_prev_HIV`, `gamma`) |
| `Value` | Calibrated or fixed value of the parameter |
| `Description` | Human-readable description of the parameter |

**Rows:** ~28 (one per parameter)

---

## Strategy extraction outputs

### `TB_COMPARTMENT_SIZES.csv` / `TB_COMPARTMENT_SIZES_OPTIMISTIC.csv`

Equilibrium compartment sizes by strategy and stratum after 50 years of simulation.

**Source:** `TB_EXTRACT_STRATEGIES.R` (both scenarios in one run; also produced as a byproduct of `extract_all_model_outputs.R` for plausible only)

| Column | Description |
|---|---|
| `Scenario` | Strategy name (or "Baseline") |
| `Stratum` | Risk stratum (HIV / Medical / FB / USB) |
| `N_pop` | Total population in that stratum |
| `S` | Susceptible compartment size |
| `L_f`, `L_s` | Fast and slow latent (unvaccinated) |
| `L_fv`, `L_sv` | Fast and slow latent (vaccinated) |
| `I` | Infectious (unvaccinated) |
| `I_v` | Infectious (vaccinated) |
| `R` | Recovered/post-treatment |
| `Total_LTBI` | Sum of all 4 latent compartments |
| `Total_LTBI_unvax` | `L_f + L_s` |
| `Total_LTBI_vax` | `L_fv + L_sv` |
| `Total_Infectious` | `I + I_v` |
| `Total_N` | Sum of all compartments |
| `LTBI_prevalence_pct` | LTBI as % of total stratum population |
| `Lf_fraction_of_LTBI` | Fast latent as % of total LTBI |
| `psi` | Vaccination rate used |
| `theta_HIV`, `theta_Medical`, `theta_FB`, `theta_USB` | Targeting weights by stratum |
| `VE` | Vaccine efficacy used |

**Rows:** 40 (8 strategies × 4 strata, plus Baseline × 4 strata in plausible; 7 strategies for optimistic)

### `TB_STRATEGY_SUMMARY.csv` / `TB_STRATEGY_SUMMARY_OPTIMISTIC.csv`

Per-strategy summary statistics at equilibrium.

**Source:** `TB_EXTRACT_STRATEGIES.R` (both scenarios in one run; also produced by `extract_all_model_outputs.R` for plausible only)

| Column | Description |
|---|---|
| `Scenario` | Strategy name |
| `psi` | Vaccination rate used |
| `VE` | Vaccine efficacy used |
| `Baseline_cases` | Annual TB cases under no-vaccination baseline |
| `Annual_cases` | Annual TB cases under this strategy |
| `Cases_prevented` | `Baseline_cases - Annual_cases` |
| `Pct_reduction` | Cases prevented as % of baseline |
| `Annual_vaccinations` | Annual vaccination doses delivered |
| `NNV` | Number needed to vaccinate to prevent one case |

**Rows:** 8 (one per strategy including Baseline)

### `TB_CASE_SOURCES.csv` / `TB_CASE_SOURCES_OPTIMISTIC.csv`

Decomposition of TB cases into transmission pathway sources at equilibrium.

**Source:** `TB_EXTRACT_STRATEGIES.R` (both scenarios; also produced by `extract_all_model_outputs.R` for plausible only)

| Column | Description |
|---|---|
| `Scenario` | Strategy name |
| `Total_cases` | Annual TB cases (all sources combined) |
| `From_Lf` | Cases arising from fast progression (recent infection) |
| `From_Ls` | Cases arising from reactivation (slow latent) |
| `From_relapse` | Cases arising from post-treatment relapse |
| `Pct_from_Lf` | `From_Lf / Total_cases × 100` |
| `Pct_from_Ls` | `From_Ls / Total_cases × 100` |
| `Pct_from_relapse` | `From_relapse / Total_cases × 100` |
| `Pct_reactivation` | `Pct_from_Ls + Pct_from_relapse` (reactivation-dominated metric) |

**Rows:** 8

### `TB_LTBI_POOL_ANALYSIS.csv` / `TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv`

Aggregated LTBI pool changes by strategy (summed across strata).

**Source:** `TB_EXTRACT_STRATEGIES.R` (both scenarios; also produced by `extract_all_model_outputs.R` for plausible only)

| Column | Description |
|---|---|
| `Scenario` | Strategy name |
| `L_f`, `L_s`, `L_fv`, `L_sv` | Aggregated fast/slow × unvax/vax compartments |
| `Total_LTBI` | Sum of all four latent compartments |
| `Total_LTBI_unvax` | `L_f + L_s` summed |
| `Total_LTBI_vax` | `L_fv + L_sv` summed |
| `Total_Infectious` | `I + I_v` summed |
| `S`, `R` | Susceptible and recovered totals |
| `LTBI_pct_of_baseline` | Total LTBI as % of baseline LTBI |
| `Lf_pct_of_baseline` | Fast latent as % of baseline fast latent |
| `Ls_pct_of_baseline` | Slow latent as % of baseline slow latent |
| `Pct_LTBI_vaccinated` | `Total_LTBI_vax / Total_LTBI × 100` |
| `Lf_fraction` | `L_f / Total_LTBI × 100` |
| `Infectious_pct_of_baseline` | Infectious as % of baseline infectious |

**Rows:** 8

### `TB_ANNUAL_FLOWS.csv` / `TB_ANNUAL_FLOWS_OPTIMISTIC.csv`

Annual transition flows between compartments at equilibrium (diagnostic / not used in figures).

**Source:** `TB_EXTRACT_STRATEGIES.R` (both scenarios; also produced by `extract_all_model_outputs.R` for plausible only)

| Column | Description |
|---|---|
| `Stratum` | Risk stratum |
| `Scenario` | Strategy name |
| `Flow_Lf_to_I`, `Flow_Ls_to_I` | New cases from unvaccinated latent compartments |
| `Flow_Lfv_to_I`, `Flow_Lsv_to_I` | New cases from vaccinated latent compartments |
| `Flow_R_to_I` | Relapse cases |
| `Flow_Lf_to_Lfv`, `Flow_Ls_to_Lsv` | New vaccinations (entering vaccinated state) |
| `Flow_Lfv_to_Lf`, `Flow_Lsv_to_Ls` | Vaccine waning (exiting vaccinated state) |
| `Flow_Lf_to_Ls`, `Flow_Lfv_to_Lsv` | Fast → slow latent transitions |
| `Flow_Lf_clear`, `Flow_Ls_clear` | Spontaneous LTBI clearance |
| `Flow_I_to_R`, `Flow_Iv_to_R` | Treatment completions |
| `Flow_R_to_Ls` | Post-treatment stabilization |
| `Total_new_cases` | Sum of all `_to_I` flows |
| `Total_vaccinations` | Sum of all `_to_Lfv` and `_to_Lsv` flows |

**Rows:** 32 (8 strategies × 4 strata)

### `TB_MODEL_OUTPUTS_COMPLETE.csv` / `TB_MODEL_OUTPUTS_COMPLETE_OPTIMISTIC.csv`

Master merge of strategy summary + LTBI pool + case sources.

**Source:** `TB_EXTRACT_MASTER.R` (post-processes outputs from `TB_EXTRACT_STRATEGIES.R`; also produced by `extract_all_model_outputs.R` for plausible only)

**Columns:** Union of `TB_STRATEGY_SUMMARY`, `TB_LTBI_POOL_ANALYSIS`, and `TB_CASE_SOURCES` columns (without duplicates), keyed on `Scenario`.

**Rows:** 8

---

## Stratum-level impact

### `TB_STRATUM_LEVEL_IMPACT.csv` / `TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv`

Per-stratum cases prevented and NNV at equilibrium for each strategy.

**Sources:**
- Plausible: `TB_STRATUM_IMPACT_PLAUSIBLE.R` or `TB_STRATUM_IMPACT.R` (also `extract_all_model_outputs.R`)
- Optimistic: `TB_STRATUM_IMPACT_OPTIMISTIC.R` or `TB_STRATUM_IMPACT.R`

| Column | Description |
|---|---|
| `Scenario` | Strategy name |
| `Stratum` | Risk stratum (HIV / Medical / FB / USB) |
| `Baseline_cases` | Annual cases in this stratum under baseline |
| `With_vaccination` | Annual cases in this stratum under this strategy |
| `Cases_prevented` | `Baseline_cases - With_vaccination` |
| `Pct_reduction` | Prevented as % of baseline |
| `Annual_vaccinations` | Annual vaccinations delivered to this stratum |
| `NNV` | Stratum-specific NNV |
| `Incidence_per_100k` | Stratum incidence rate per 100,000 under strategy |
| `Baseline_incidence` | Baseline incidence rate per 100,000 |

**Rows:** 40 (plausible: 9 strategies + Baseline × 4 strata) or 36 (optimistic: 8 strategies + Baseline × 4 strata)

---

## Time-to-impact

### `TB_TIME_TO_IMPACT_RESULTS.csv` / `TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv`

Annual and cumulative impact at years 1, 2, 3, 5, 10, 15, 20, 25, 30.

**Sources:**
- Plausible: `TB_TIME_TO_IMPACT_PLAUSIBLE.R` or `TB_TIME_TO_IMPACT.R`
- Optimistic: `TB_TIME_TO_IMPACT_OPTIMISTIC.R` or `TB_TIME_TO_IMPACT.R`

| Column | Description |
|---|---|
| `year` | Year of follow-up (1-30) |
| `strategy` | Strategy name |
| `baseline_cases` | Annual baseline cases at this year (scaled to 2024) |
| `intervention_cases` | Annual cases under intervention at this year |
| `cases_prevented` | `baseline_cases - intervention_cases` (annual rate) |
| `pct_reduction` | Annual % reduction at this year |
| `annual_vax` | Annual vaccinations delivered at this year |
| `cumul_prevented` | Cumulative cases prevented from year 0 to this year |
| `cumul_vax` | Cumulative vaccinations from year 0 to this year |
| `cumul_NNV` | `cumul_vax / cumul_prevented` |
| `prev_HIV`, `prev_Med`, `prev_FB`, `prev_USB` | Cases prevented per year by stratum |

**Rows:** 63 (7 strategies × 9 time points) for optimistic; 54 (6 × 9) for plausible

---

## Effective coverage

### `TB_EFFECTIVE_COVERAGE_OVERALL.csv` / `_OPTIMISTIC.csv`

Effective vaccine coverage across the entire infected population at each time point.

**Sources:**
- Plausible: `TB_EFFECTIVE_COVERAGE_PLAUSIBLE.R` or `TB_EFFECTIVE_COVERAGE.R`
- Optimistic: `TB_EFFECTIVE_COVERAGE_OPTIMISTIC.R` or `TB_EFFECTIVE_COVERAGE.R`

| Column | Description |
|---|---|
| `year` | Year (1, 2, 3, 5, 10, 15, 20, 25, 30) |
| `strategy` | Strategy name |
| `eff_coverage_pct` | `(L_fv + L_sv) / (L_f + L_s + L_fv + L_sv) × 100`, summed across all 4 strata |

**Rows:** 54 (6 strategies × 9 time points)

### `TB_EFFECTIVE_COVERAGE_STRATUM.csv` / `_OPTIMISTIC.csv`

Effective coverage broken down by stratum.

**Sources:** Same as `TB_EFFECTIVE_COVERAGE_OVERALL.csv` (above).

| Column | Description |
|---|---|
| `year` | Year |
| `strategy` | Strategy name |
| `stratum` | Risk stratum |
| `eff_coverage_pct` | Coverage within that stratum |
| `targeted` | TRUE/FALSE — whether this stratum is targeted by the strategy |
| `eff_coverage_plot` | Coverage with NA for non-targeted (used for plotting masks) |

**Rows:** 216 (6 × 4 × 9) or 252 (7 × 4 × 9 for older optimistic version)

### `TB_EFFECTIVE_COVERAGE_TARGETED.csv` / `_OPTIMISTIC.csv`

Effective coverage pooled across only the strategy's targeted strata.

**Sources:** Same as `TB_EFFECTIVE_COVERAGE_OVERALL.csv` (above).

| Column | Description |
|---|---|
| `year` | Year |
| `strategy` | Strategy name |
| `eff_coverage_pct` | Coverage within targeted strata combined |

**Rows:** 54

### `TB_EFFECTIVE_COVERAGE_POOLED_OPTIMISTIC.csv` (legacy)

Equilibrium-only summary table (pre-modular).

**Source:** Legacy file from earlier optimistic analysis work (no script in current pipeline regenerates it; superseded by `TB_EFFECTIVE_COVERAGE.R` outputs).

| Column | Description |
|---|---|
| `Strategy` | Strategy name |
| `Targeted_pct` | Equilibrium effective coverage in targeted strata |
| `Overall_pct` | Equilibrium effective coverage across all strata |

**Rows:** 6

---

## Direct/indirect + psi sensitivity + threshold

### `TB_DIRECT_INDIRECT_EFFECTS.csv` / `_OPTIMISTIC.csv`

Decomposition of total cases prevented into direct (vaccinated person protected) vs. indirect (transmission reduction) effects, computed via static-FOI counterfactual.

**Sources:**
- Plausible: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_PLAUSIBLE.R` or `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`
- Optimistic: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_OPTIMISTIC.R`, `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`, or `TB_DIRECT_INDIRECT_OPTIMISTIC.R` (legacy subset)

| Column | Description |
|---|---|
| `Strategy` | Strategy name |
| `Total_prevented` | Total cases prevented under full dynamic model |
| `Direct_prevented` | Cases prevented under static-FOI counterfactual (direct effect only) |
| `Indirect_prevented` | `Total_prevented - Direct_prevented` |
| `Pct_direct` | `Direct_prevented / Total_prevented × 100` |
| `Pct_indirect` | `Indirect_prevented / Total_prevented × 100` |
| `Total_pct_reduction` | `Total_prevented / baseline_cases × 100` |

**Rows:** 8 (plausible) or 7 (optimistic)

### `TB_PSI_SENSITIVITY.csv` / `_OPTIMISTIC.csv`

Vaccination rate (psi) dose-response sweep.

**Sources:**
- Plausible: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_PLAUSIBLE.R` or `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`
- Optimistic: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_OPTIMISTIC.R` or `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`

| Column | Description |
|---|---|
| `strategy` | Strategy name |
| `psi_pct` | psi value as % per year (1, 2, 3, 5, 7, 10, 15, 20) |
| `annual_cases` | Annual cases at this psi |
| `cases_prevented` | Cases prevented vs. baseline |
| `pct_reduction` | % reduction |
| `annual_vax` | Annual vaccinations |
| `NNV` | Standard NNV |
| `marginal_NNV` | Vaccinations to prevent one additional case at this psi step (marginal efficiency) |

**Rows:** 64 (8 strategies × 8 psi values)

### `TB_THRESHOLD_ANALYSIS.csv` / `_OPTIMISTIC.csv`

Full grid of vaccine efficacy × duration × strategy outcomes.

**Sources:**
- Plausible: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_PLAUSIBLE.R` or `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`
- Optimistic: `TB_DIRECT_INDIRECT_PSI_THRESHOLD_OPTIMISTIC.R` or `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R`

| Column | Description |
|---|---|
| `strategy` | Strategy name |
| `VE` | Vaccine efficacy as % (10-95 in steps of 5) |
| `duration_yr` | Duration of protection in years (5, 10, 15, 20, 30) |
| `cases_prevented` | Cases prevented at this combination |
| `pct_reduction` | % reduction |
| `NNV` | NNV at this combination |

**Rows:** 720 (8 strategies × 18 VE × 5 durations)

### `TB_THRESHOLD_MINVE.csv` / `_OPTIMISTIC.csv`

Minimum VE required to achieve specific clinical targets at each duration.

**Sources:** Same as `TB_THRESHOLD_ANALYSIS.csv` (above).

| Column | Description |
|---|---|
| `strategy` | Strategy name |
| `duration_yr` | Duration in years |
| `target` | Target description (e.g., "≥5% reduction", "NNV < 500") |
| `min_VE_pct` | Minimum VE (%) to meet the target, or NA if unattainable |

**Rows:** 160 (8 strategies × 5 durations × 4 targets)

---

## One-way sensitivity + published NNV

### `TB_ONEWAY_SENSITIVITY.csv` / `_OPTIMISTIC.csv`

Per-parameter sweep with recalibration at each value.

**Sources:**
- Plausible: `TB_ONEWAY_PLAUSIBLE.R` or `TB_ONEWAY.R`
- Optimistic: `TB_ONEWAY_OPTIMISTIC.R` or `TB_ONEWAY.R`

| Column | Description |
|---|---|
| `parameter` | Parameter being varied (e.g., `gamma`, `eps_s_USB`) |
| `value` | Parameter value at this draw |
| `is_baseline` | TRUE if `value` matches the baseline value |
| `cal_beta` | Recalibrated transmission rate β at this draw |
| `cases_prevented` | Cases prevented at All High-Risk reference strategy |
| `pct_reduction` | % reduction |
| `NNV` | NNV at this parameter value |
| `baseline_cases` | 2024 actual cases (constant reference) |

**Rows:** 64 (varies by parameter; ~6-8 values per parameter × 11 parameters)

### `TB_ONEWAY_TORNADO.csv` / `_OPTIMISTIC.csv`

Range summary per parameter (used for tornado plots).

**Sources:** Same as `TB_ONEWAY_SENSITIVITY.csv` (above).

| Column | Description |
|---|---|
| `parameter` | Parameter name |
| `label` | Human-readable label |
| `baseline_prevented` | Cases prevented at baseline parameter value |
| `prevented_at_lo` | Cases prevented at low end of parameter range |
| `prevented_at_hi` | Cases prevented at high end of parameter range |
| `range_prevented` | `max - min` of cases prevented across the range |
| `lo_value` | Low end of parameter range |
| `hi_value` | High end of parameter range |

**Rows:** 11-12 (one per parameter swept)

### `TB_PUBLISHED_NNV_COMPARISON.csv`

Literature comparison of NNV across published TB vaccine modeling studies. Plausible-only (literature values are scenario-independent).

**Source:** `TB_ONEWAY_PLAUSIBLE.R` or `TB_ONEWAY.R`

| Column | Description |
|---|---|
| `Study` | Citation (author and year) |
| `Setting` | Country/region of analysis |
| `Population` | Target population (e.g., "Adolescent/adult, post-infection") |
| `VE_assumed` | VE assumption used in the study |
| `Duration` | Duration of protection assumed |
| `NNV_reported` | Published NNV range or point estimate |
| `Notes` | Methodological context |

**Rows:** 8

---

## LHS uncertainty + PRCC

### `TB_LHS_UNCERTAINTY_INTERVALS.csv` / `_OPTIMISTIC.csv`

95% uncertainty intervals from 500 LHS draws.

**Sources:**
- Plausible: `TB_LHS_PLAUSIBLE.R` or `TB_LHS.R`
- Optimistic: `TB_LHS_OPTIMISTIC.R` or `TB_LHS.R`

| Column | Description |
|---|---|
| `Strategy` | Strategy name |
| `Prevented_median` | Median cases prevented across valid draws |
| `Prevented_lo`, `Prevented_hi` | 2.5th and 97.5th percentiles |
| `PctRed_median` | Median % reduction |
| `PctRed_lo`, `PctRed_hi` | 2.5th and 97.5th percentiles |
| `NNV_median` | Median NNV |
| `NNV_lo`, `NNV_hi` | 2.5th and 97.5th percentiles |
| `N_valid` | Number of LHS draws (out of 500) where calibration succeeded and cases_prevented > 0 |

**Rows:** 6-7 (one per primary strategy)

### `TB_PRCC_RESULTS.csv` / `_OPTIMISTIC.csv`

Partial rank correlation coefficients for each parameter against outcome metrics.

**Sources:** Same as `TB_LHS_UNCERTAINTY_INTERVALS.csv` (above).

| Column | Description |
|---|---|
| `parameter` | Parameter name |
| `PRCC` | Partial rank correlation coefficient (-1 to +1) |
| `p_value` | Significance test p-value |
| `outcome` | Outcome metric (`cases_prevented` or `NNV`) |

**Rows:** 24 (12 parameters × 2 outcomes for plausible) or 12 (12 × 1 outcome for optimistic, cases_prevented only)

---

## NNV benchmarking

### `TB_NNV_BENCHMARK_OPTIMISTIC.csv`

Equilibrium NNV point estimates for the cross-vaccine benchmark figure (optimistic only — plausible NNV values come from `TB_STRATEGY_SUMMARY.csv`).

**Source:** `TB_NNV_BENCHMARK_OPTIMISTIC.R`

| Column | Description |
|---|---|
| `Strategy` | Strategy name |
| `VE` | Vaccine efficacy used (0.70) |
| `psi` | Vaccination rate used (0.50) |
| `Cases_prevented` | Annual cases prevented at equilibrium |
| `Annual_vaccinations` | Annual vaccinations delivered |
| `NNV` | Number needed to vaccinate |
| `Scenario` | Always "Optimistic" |

**Rows:** 7

---

## Alternative calibration

### `TB_CLEAN_ALTERNATIVE_CALIBRATION_PARAMS.csv`

Comparison of primary vs. alternative calibrated parameter values.

**Source:** `TB_CLEAN_ALTERNATIVE_CALIBRATION.R` (project folder, not in outputs)

| Column | Description |
|---|---|
| `Parameter` | Parameter name |
| `Primary` | Value under primary calibration (anchored on Ekramnia reactivation rates) |
| `Alternative` | Value under alternative calibration (anchored on Aim 1 ecological prevalence) |
| `Ratio` | `Alternative / Primary` |

**Rows:** 11

### `TB_CLEAN_ALTERNATIVE_INTERVENTION_RESULTS.csv`

Per-strategy intervention outcomes under both calibrations.

**Source:** `TB_CLEAN_ALTERNATIVE_CALIBRATION.R` (project folder, not in outputs)

| Column | Description |
|---|---|
| `strategy` | Strategy name |
| `calibration` | "primary" or "alternative" |
| `baseline_cases` | Baseline cases under this calibration |
| `intervention_cases` | Cases under this strategy |
| `cases_prevented` | Cases prevented |
| `pct_reduction` | % reduction |
| `annual_vaccinations` | Annual vaccinations |
| `NNV` | NNV |

**Rows:** 14 (7 strategies × 2 calibrations)

### `TB_CLEAN_ALTERNATIVE_SUMMARY.csv`

Side-by-side comparison of primary and alternative calibration outcomes.

**Source:** `TB_CLEAN_ALTERNATIVE_CALIBRATION.R` (project folder, not in outputs)

| Column | Description |
|---|---|
| `Strategy` | Strategy name |
| `Primary_Prevented`, `Primary_NNV`, `Primary_PctRed` | Outcomes under primary calibration |
| `Alt_Prevented`, `Alt_NNV`, `Alt_PctRed` | Outcomes under alternative calibration |
| `Prevented_Ratio` | `Alt_Prevented / Primary_Prevented` (sanity check; ~1.0 means similar) |
| `NNV_Ratio` | `Alt_NNV / Primary_NNV` |

**Rows:** 7

---

## Legacy / supporting files

### `TB_OPTIMISTIC_SCENARIOS.csv`

Multi-psi/multi-VE sweep used for the appendix coverage analysis. Has effective coverage columns by stratum at each parameter combination.

**Source:** Legacy file produced by an earlier version of the optimistic appendix-extension scripts; not regenerated by the current modular pipeline. Consumed by Figure 12 (threshold heatmap) and Figure 14 (psi dose-response).

| Column | Description |
|---|---|
| `strategy` | Strategy name |
| `psi_pct` | psi value as % (e.g., 1, 2, 5, 10, 20, 50) |
| `VE_pct` | VE as % (e.g., 50, 70) |
| `annual_cases` | Annual cases |
| `cases_prevented` | Cases prevented |
| `pct_reduction` | % reduction |
| `annual_vax` | Annual vaccinations |
| `NNV` | NNV |
| `eff_cov_overall` | Equilibrium effective coverage across all strata |
| `eff_cov_targeted` | Equilibrium effective coverage in targeted strata |
| `eff_cov_HIV`, `eff_cov_Medical`, `eff_cov_NUSB`, `eff_cov_USB` | Stratum-specific effective coverage |

**Rows:** 90 (15 psi × 6 strategies)

### `TB_EXPANDED_COVERAGE_ANALYSIS.csv`

Same column structure as `TB_OPTIMISTIC_SCENARIOS.csv` but with an added `duration_yr` column for duration sensitivity. **Rows:** 210.

**Source:** Legacy file from earlier appendix coverage analysis; not regenerated by the current modular pipeline.

### `TB_LONG_DURATION_SCENARIOS.csv`

Same as `TB_EXPANDED_COVERAGE_ANALYSIS.csv` for long-duration variants. **Rows:** 48.

**Source:** Legacy file from earlier long-duration sensitivity work; not regenerated by the current modular pipeline.

### `TB_PSI_SENSITIVITY_EXTENDED.csv`

Extended psi sweep at VE=50% (12 psi values from 1% to 50%) with effective coverage by stratum. Same column structure as `TB_OPTIMISTIC_SCENARIOS.csv`. **Rows:** 72 (12 psi × 6 strategies).

**Source:** Legacy file from extended psi sweep work; not regenerated by the current modular pipeline. The current `TB_PSI_SENSITIVITY.csv` covers the standard 1-20% psi range used in the manuscript.

---

## Quick reference: which figures consume which CSVs?

The manuscript figures are produced by `TB_DISSERTATION_FIGURES_4_22_26_v10.R`. Many CSVs are loaded once at the top of the script and reused across multiple figures; others are loaded inline within a specific figure block. The table below maps each script-internal figure number to its manuscript figure label and the CSVs it consumes.

| Script Fig # | Manuscript label | Description | CSVs consumed |
|---|---|---|---|
| 1 | (model diagram) | Model structure diagram | None — drawn programmatically |
| 2 | Figure S1 | Case source validation | `TB_CASE_SOURCES.csv`, `TB_ANNUAL_FLOWS.csv` |
| 3 | **Figure 2** | Strategy comparison (cases prevented + NNV, both scenarios) | `TB_STRATEGY_SUMMARY.csv`, `TB_OPTIMISTIC_SCENARIOS.csv`, `TB_THRESHOLD_ANALYSIS.csv` |
| 4 | Figure S2 | Stratum-level impact breakdown | `TB_STRATUM_LEVEL_IMPACT.csv`, `TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv` |
| 5 | (auxiliary) | Incidence trajectory over time | `TB_TIME_TO_IMPACT_RESULTS.csv` |
| 6 | (auxiliary) | Cumulative cases prevented over time | `TB_TIME_TO_IMPACT_RESULTS.csv` |
| 7 | (auxiliary) | Time-to-impact curves (annual) | `TB_TIME_TO_IMPACT_RESULTS.csv` |
| 8 | Figure S3 | Tornado diagram (one-way sensitivity) | `TB_ONEWAY_TORNADO.csv` |
| 9 | (auxiliary) | LHS forest plot | `TB_LHS_UNCERTAINTY_INTERVALS.csv` |
| 10 | Figure S4 | PRCC bar plot (dual-scenario) | `TB_PRCC_RESULTS.csv`, `TB_PRCC_RESULTS_OPTIMISTIC.csv` |
| 11 | (auxiliary) | Direct vs. indirect effects | `TB_DIRECT_INDIRECT_EFFECTS.csv` |
| 12 | **Figure 3** | Threshold heatmap (VE × duration, dual-panel) | `TB_THRESHOLD_ANALYSIS.csv`, `TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv` |
| 13 | Figure S5 | NNV comparison with established vaccines | `TB_STRATEGY_SUMMARY.csv`, `TB_LHS_UNCERTAINTY_INTERVALS.csv`, `TB_OPTIMISTIC_SCENARIOS.csv`, `TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv` |
| 14 | (auxiliary) | Psi dose-response | `TB_PSI_SENSITIVITY.csv` |

### Notes on this mapping

**`TB_NNV_BENCHMARK_OPTIMISTIC.csv` is not consumed by any figure.** Despite being part of the production pipeline, the figures script reads optimistic NNV point estimates from `TB_OPTIMISTIC_SCENARIOS.csv` (variable `opt_strat`) rather than the dedicated benchmark file. Keep the script in case you want a clean, focused output for downstream reporting, but the figure does not depend on it.

**Script figure numbering does not match manuscript figure numbering.** This is because the script generates many auxiliary diagnostic figures that don't appear in the manuscript or appendix. Use the "Manuscript label" column to find the figure as numbered in the paper.

**Tables S5–S7 in the appendix come from CSVs not figures.** The manuscript and appendix include several supplementary tables that report numbers from CSVs directly without a figure rendering:

| Manuscript table | CSV source |
|---|---|
| Table 2 (main) | `TB_STRATEGY_SUMMARY.csv` + `TB_LHS_UNCERTAINTY_INTERVALS.csv` + `TB_OPTIMISTIC_SCENARIOS.csv` + `TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv` |
| Table S1 | `TB_CALIBRATION_RESULTS.csv` |
| Table S4 | `TB_DIRECT_INDIRECT_EFFECTS.csv` + `TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv` |
| Table S5 | `TB_PSI_SENSITIVITY.csv` (or `TB_OPTIMISTIC_SCENARIOS.csv` for extended psi) |
| Table S6 | `TB_OPTIMISTIC_SCENARIOS.csv` (effective coverage at equilibrium) |
| Table S7 / S8 | `TB_TIME_TO_IMPACT_RESULTS.csv` + `TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv` |
