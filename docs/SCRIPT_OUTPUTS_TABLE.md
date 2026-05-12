# Script-to-Output Reference Table

This table lists every R script in the analysis pipeline and the exact files each one produces. CSVs are listed first within each row, then RData files. All output paths are relative to the working directory.

---

## Master pipeline runner

| Script | Outputs produced | Notes |
|---|---|---|
| `RUN_ALL.R` | (none directly — sources other scripts) | Executes 11 analysis scripts in sequence |

---

## Base model (project folder)

| Script | Outputs produced | Notes |
|---|---|---|
| `TB_VACCINATION_MODEL_COMPLETE.R` | (none — produces in-memory `calibration` object) | Sourced automatically by every analysis script |

---

## Output extraction

| Script | Outputs produced |
|---|---|
| `TB_EXTRACT_CALIBRATION.R` | `TB_CALIBRATION_RESULTS.csv` |
| `TB_EXTRACT_STRATEGIES.R` | `TB_COMPARTMENT_SIZES.csv` <br> `TB_STRATEGY_SUMMARY.csv` <br> `TB_CASE_SOURCES.csv` <br> `TB_LTBI_POOL_ANALYSIS.csv` <br> `TB_ANNUAL_FLOWS.csv` <br> `TB_COMPARTMENT_SIZES_OPTIMISTIC.csv` <br> `TB_STRATEGY_SUMMARY_OPTIMISTIC.csv` <br> `TB_CASE_SOURCES_OPTIMISTIC.csv` <br> `TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv` <br> `TB_ANNUAL_FLOWS_OPTIMISTIC.csv` |
| `TB_EXTRACT_MASTER.R` | `TB_MODEL_OUTPUTS_COMPLETE.csv` <br> `TB_MODEL_OUTPUTS_COMPLETE_OPTIMISTIC.csv` |

---

## Stratum-level impact

| Script | Outputs produced |
|---|---|
| `TB_STRATUM_IMPACT_PLAUSIBLE.R` | `TB_STRATUM_LEVEL_IMPACT.csv` |
| `TB_STRATUM_IMPACT_OPTIMISTIC.R` | `TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv` |
| `TB_STRATUM_IMPACT.R` | `TB_STRATUM_LEVEL_IMPACT.csv` <br> `TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv` |

---

## Time-to-impact / transient dynamics

| Script | Outputs produced |
|---|---|
| `TB_TIME_TO_IMPACT_PLAUSIBLE.R` | `TB_TIME_TO_IMPACT_RESULTS.csv` <br> `TB_TIME_TO_IMPACT_RESULTS.RData` |
| `TB_TIME_TO_IMPACT_OPTIMISTIC.R` | `TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv` |
| `TB_TIME_TO_IMPACT.R` | `TB_TIME_TO_IMPACT_RESULTS.csv` <br> `TB_TIME_TO_IMPACT_RESULTS.RData` <br> `TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv` <br> `TB_TIME_TO_IMPACT_OPTIMISTIC_RESULTS.RData` |

---

## Effective coverage

| Script | Outputs produced |
|---|---|
| `TB_EFFECTIVE_COVERAGE_PLAUSIBLE.R` | `TB_EFFECTIVE_COVERAGE_OVERALL.csv` <br> `TB_EFFECTIVE_COVERAGE_STRATUM.csv` <br> `TB_EFFECTIVE_COVERAGE_TARGETED.csv` |
| `TB_EFFECTIVE_COVERAGE_OPTIMISTIC.R` | `TB_EFFECTIVE_COVERAGE_OVERALL_OPTIMISTIC.csv` <br> `TB_EFFECTIVE_COVERAGE_STRATUM_OPTIMISTIC.csv` <br> `TB_EFFECTIVE_COVERAGE_TARGETED_OPTIMISTIC.csv` |
| `TB_EFFECTIVE_COVERAGE.R` | `TB_EFFECTIVE_COVERAGE_OVERALL.csv` <br> `TB_EFFECTIVE_COVERAGE_STRATUM.csv` <br> `TB_EFFECTIVE_COVERAGE_TARGETED.csv` <br> `TB_EFFECTIVE_COVERAGE_OVERALL_OPTIMISTIC.csv` <br> `TB_EFFECTIVE_COVERAGE_STRATUM_OPTIMISTIC.csv` <br> `TB_EFFECTIVE_COVERAGE_TARGETED_OPTIMISTIC.csv` |

---

## Direct/indirect + psi sensitivity + threshold (3-in-1)

| Script | Outputs produced |
|---|---|
| `TB_DIRECT_INDIRECT_PSI_THRESHOLD_PLAUSIBLE.R` | `TB_DIRECT_INDIRECT_EFFECTS.csv` <br> `TB_PSI_SENSITIVITY.csv` <br> `TB_THRESHOLD_ANALYSIS.csv` <br> `TB_THRESHOLD_MINVE.csv` <br> `TB_DIPT_PLAUSIBLE.RData` |
| `TB_DIRECT_INDIRECT_PSI_THRESHOLD_OPTIMISTIC.R` | `TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv` <br> `TB_PSI_SENSITIVITY_OPTIMISTIC.csv` <br> `TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv` <br> `TB_THRESHOLD_MINVE_OPTIMISTIC.csv` <br> `TB_DIPT_OPTIMISTIC.RData` |
| `TB_DIRECT_INDIRECT_PSI_THRESHOLD.R` | `TB_DIRECT_INDIRECT_EFFECTS.csv` <br> `TB_PSI_SENSITIVITY.csv` <br> `TB_THRESHOLD_ANALYSIS.csv` <br> `TB_THRESHOLD_MINVE.csv` <br> `TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv` <br> `TB_PSI_SENSITIVITY_OPTIMISTIC.csv` <br> `TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv` <br> `TB_THRESHOLD_MINVE_OPTIMISTIC.csv` <br> `TB_DIPT.RData` |
| `TB_DIRECT_INDIRECT_OPTIMISTIC.R` (subset utility) | `TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv` |

---

## One-way sensitivity + published NNV

| Script | Outputs produced |
|---|---|
| `TB_ONEWAY_PLAUSIBLE.R` | `TB_ONEWAY_SENSITIVITY.csv` <br> `TB_ONEWAY_TORNADO.csv` <br> `TB_PUBLISHED_NNV_COMPARISON.csv` <br> `TB_ONEWAY_PLAUSIBLE.RData` |
| `TB_ONEWAY_OPTIMISTIC.R` | `TB_ONEWAY_SENSITIVITY_OPTIMISTIC.csv` <br> `TB_ONEWAY_TORNADO_OPTIMISTIC.csv` <br> `TB_ONEWAY_OPTIMISTIC.RData` |
| `TB_ONEWAY.R` | `TB_ONEWAY_SENSITIVITY.csv` <br> `TB_ONEWAY_TORNADO.csv` <br> `TB_PUBLISHED_NNV_COMPARISON.csv` <br> `TB_ONEWAY_SENSITIVITY_OPTIMISTIC.csv` <br> `TB_ONEWAY_TORNADO_OPTIMISTIC.csv` <br> `TB_ONEWAY.RData` |

---

## LHS uncertainty + PRCC

| Script | Outputs produced |
|---|---|
| `TB_LHS_PLAUSIBLE.R` | `TB_LHS_UNCERTAINTY_INTERVALS.csv` <br> `TB_PRCC_RESULTS.csv` <br> `TB_LHS_RESULTS.RData` |
| `TB_LHS_OPTIMISTIC.R` | `TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv` <br> `TB_PRCC_RESULTS_OPTIMISTIC.csv` <br> `TB_LHS_OPTIMISTIC_RESULTS.RData` |
| `TB_LHS.R` | `TB_LHS_UNCERTAINTY_INTERVALS.csv` <br> `TB_PRCC_RESULTS.csv` <br> `TB_LHS_RESULTS.RData` <br> `TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv` <br> `TB_PRCC_RESULTS_OPTIMISTIC.csv` <br> `TB_LHS_OPTIMISTIC_RESULTS.RData` |

---

## NNV benchmarking

| Script | Outputs produced |
|---|---|
| `TB_NNV_BENCHMARK_OPTIMISTIC.R` | `TB_NNV_BENCHMARK_OPTIMISTIC.csv` |

(Plausible NNV values come from `TB_STRATEGY_SUMMARY.csv` and `TB_LHS_UNCERTAINTY_INTERVALS.csv`, both produced earlier in the pipeline. No dedicated plausible script.)

---

## Alternative calibration (project folder)

| Script | Outputs produced |
|---|---|
| `TB_CLEAN_ALTERNATIVE_CALIBRATION.R` | `TB_CLEAN_ALTERNATIVE_CALIBRATION_PARAMS.csv` <br> `TB_CLEAN_ALTERNATIVE_INTERVENTION_RESULTS.csv` <br> `TB_CLEAN_ALTERNATIVE_SUMMARY.csv` |

---

## Full pipeline output count

Running `RUN_ALL.R` end-to-end (which uses unified scripts where available) produces approximately 50 unique output files: roughly 45 CSVs and 5 RData files. Listed by scenario:

**Scenario-independent (1 file)** <br>
`TB_CALIBRATION_RESULTS.csv`

**Plausible scenario (~22 files)** <br>
Strategy outputs: `TB_COMPARTMENT_SIZES.csv`, `TB_STRATEGY_SUMMARY.csv`, `TB_CASE_SOURCES.csv`, `TB_LTBI_POOL_ANALYSIS.csv`, `TB_ANNUAL_FLOWS.csv`, `TB_MODEL_OUTPUTS_COMPLETE.csv`, `TB_STRATUM_LEVEL_IMPACT.csv` <br>
Time-to-impact: `TB_TIME_TO_IMPACT_RESULTS.csv` + `.RData` <br>
Effective coverage: `TB_EFFECTIVE_COVERAGE_OVERALL.csv`, `TB_EFFECTIVE_COVERAGE_STRATUM.csv`, `TB_EFFECTIVE_COVERAGE_TARGETED.csv` <br>
Direct/indirect/psi/threshold: `TB_DIRECT_INDIRECT_EFFECTS.csv`, `TB_PSI_SENSITIVITY.csv`, `TB_THRESHOLD_ANALYSIS.csv`, `TB_THRESHOLD_MINVE.csv` <br>
One-way: `TB_ONEWAY_SENSITIVITY.csv`, `TB_ONEWAY_TORNADO.csv`, `TB_PUBLISHED_NNV_COMPARISON.csv` <br>
LHS: `TB_LHS_UNCERTAINTY_INTERVALS.csv`, `TB_PRCC_RESULTS.csv` + `TB_LHS_RESULTS.RData`

**Optimistic scenario (~21 files)** <br>
Strategy outputs: `TB_COMPARTMENT_SIZES_OPTIMISTIC.csv`, `TB_STRATEGY_SUMMARY_OPTIMISTIC.csv`, `TB_CASE_SOURCES_OPTIMISTIC.csv`, `TB_LTBI_POOL_ANALYSIS_OPTIMISTIC.csv`, `TB_ANNUAL_FLOWS_OPTIMISTIC.csv`, `TB_MODEL_OUTPUTS_COMPLETE_OPTIMISTIC.csv`, `TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv` <br>
Time-to-impact: `TB_TIME_TO_IMPACT_RESULTS_OPTIMISTIC.csv` + `.RData` <br>
Effective coverage: `TB_EFFECTIVE_COVERAGE_OVERALL_OPTIMISTIC.csv`, `TB_EFFECTIVE_COVERAGE_STRATUM_OPTIMISTIC.csv`, `TB_EFFECTIVE_COVERAGE_TARGETED_OPTIMISTIC.csv` <br>
Direct/indirect/psi/threshold: `TB_DIRECT_INDIRECT_EFFECTS_OPTIMISTIC.csv`, `TB_PSI_SENSITIVITY_OPTIMISTIC.csv`, `TB_THRESHOLD_ANALYSIS_OPTIMISTIC.csv`, `TB_THRESHOLD_MINVE_OPTIMISTIC.csv` <br>
One-way: `TB_ONEWAY_SENSITIVITY_OPTIMISTIC.csv`, `TB_ONEWAY_TORNADO_OPTIMISTIC.csv` <br>
LHS: `TB_LHS_UNCERTAINTY_INTERVALS_OPTIMISTIC.csv`, `TB_PRCC_RESULTS_OPTIMISTIC.csv` + `TB_LHS_OPTIMISTIC_RESULTS.RData` <br>
NNV benchmark: `TB_NNV_BENCHMARK_OPTIMISTIC.csv`

**Combined RData files (1 file)** <br>
`TB_DIPT.RData` (unified direct/indirect/psi/threshold), `TB_ONEWAY.RData` (unified one-way)

**Alternative calibration (3 files)** <br>
`TB_CLEAN_ALTERNATIVE_CALIBRATION_PARAMS.csv`, `TB_CLEAN_ALTERNATIVE_INTERVENTION_RESULTS.csv`, `TB_CLEAN_ALTERNATIVE_SUMMARY.csv`
