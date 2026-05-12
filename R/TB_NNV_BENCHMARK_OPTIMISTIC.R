################################################################################
#
#  TB VACCINATION MODEL — NNV BENCHMARK POINT ESTIMATES (OPTIMISTIC SCENARIO)
#
################################################################################
#
# PURPOSE:
#   Compute equilibrium NNV (number needed to vaccinate to prevent one case)
#   point estimates for each strategy under the optimistic scenario. Used as
#   the optimistic-side data points in the cross-vaccine NNV benchmark figure
#   (Figure S5).
#
#   Optimistic scenario: VE = 0.70, psi = 0.50/yr for all strategies,
#   10-year duration of protection.
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_NNV_BENCHMARK_OPTIMISTIC.R")
#   3. Result saved as CSV in working directory
#
# OUTPUT FILES:
#   - TB_NNV_BENCHMARK_OPTIMISTIC.csv
#
# NOTE ON STRATEGY SET:
#   Uses 7 strategies matching the manuscript's primary + exploratory set
#   under the optimistic scenario. The "PLWH + Non-U.S.-Born" strategy
#   (theta = c(1, 0.18, 1, 0)) replaces "All Foreign-Born".
#
# EXPECTED RUNTIME: ~1-2 minutes
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB VACCINATION MODEL — NNV BENCHMARK (OPTIMISTIC)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n", calibration$beta))
cat(sprintf("  Actual 2024 cases = %d\n\n", actual_2024_cases))

baseline_cal <- calibration

# ---- Optimistic parameters -----------------------------------------------
OPT_VE       <- 0.70
OPT_PSI      <- 0.50
sim_years    <- 30

# ---- Strategies (raw R-code labels matching existing CSVs) ----------------
strategies_opt <- list(
  all_high_risk = list(name = "All High-Risk (HIV+Med+FB)",
                       psi = OPT_PSI, theta = c(1, 1, 1, 0)),
  plwh_nusb     = list(name = "PLWH + Non-U.S.-Born",
                       psi = OPT_PSI, theta = c(1, 0.18, 1, 0)),
  fb_stratum    = list(name = "FB Stratum Only",
                       psi = OPT_PSI, theta = c(0, 0, 1, 0)),
  universal     = list(name = "Universal",
                       psi = OPT_PSI, theta = c(1, 1, 1, 1)),
  medical       = list(name = "Medical Only",
                       psi = OPT_PSI, theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",
                       psi = OPT_PSI, theta = c(1, 0, 0, 0)),
  hiv_medical   = list(name = "HIV + Medical",
                       psi = OPT_PSI, theta = c(1, 1, 0, 0))
)

# ============================================================================
# RUN NNV BENCHMARK
# ============================================================================
cat("================================================================\n")
cat("  RUNNING NNV BENCHMARK (OPTIMISTIC)\n")
cat("================================================================\n\n")

# Baseline (no vaccination) for scaling factor
params_base_opt <- list(
  beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
  psi = 0, theta = c(0, 0, 0, 0),
  mtb_prev = baseline_cal$mtb_prev,
  M = M_default, N = N_vec, VE = OPT_VE
)
base_run_opt <- run_model(params_base_opt, years = sim_years)
base_cases_raw <- calc_annual_cases(base_run_opt$final_state, params_base_opt)$total
sf <- actual_2024_cases / base_cases_raw

nnv_bench_opt <- data.frame(
  Strategy = character(),
  VE = numeric(),
  psi = numeric(),
  Cases_prevented = numeric(),
  Annual_vaccinations = numeric(),
  NNV = numeric(),
  Scenario = character(),
  stringsAsFactors = FALSE
)

for (s in names(strategies_opt)) {
  strat <- strategies_opt[[s]]
  cat(sprintf("  %-30s ... ", strat$name))

  params_int <- list(
    beta = baseline_cal$beta, kappa = c(1, 1, 1, 1),
    psi = strat$psi, theta = strat$theta,
    mtb_prev = baseline_cal$mtb_prev,
    M = M_default, N = N_vec, VE = OPT_VE
  )
  run <- run_model(params_int, years = sim_years)
  cases <- calc_annual_cases(run$final_state, params_int)$total * sf
  vax   <- calc_vaccinations(run$final_state, params_int$psi, params_int$theta)$total
  prev  <- actual_2024_cases - cases
  nnv   <- if (prev > 0) vax / prev else NA

  nnv_bench_opt <- rbind(nnv_bench_opt, data.frame(
    Strategy = strat$name,
    VE = OPT_VE,
    psi = OPT_PSI,
    Cases_prevented = round(prev),
    Annual_vaccinations = round(vax),
    NNV = round(nnv),
    Scenario = "Optimistic",
    stringsAsFactors = FALSE
  ))
  cat(sprintf("NNV=%d, %d cases prevented\n", round(nnv), round(prev)))
}

write.csv(nnv_bench_opt, "TB_NNV_BENCHMARK_OPTIMISTIC.csv", row.names = FALSE)

cat("\n================================================================\n")
cat("  NNV BENCHMARK COMPLETE (OPTIMISTIC)\n")
cat("================================================================\n")
cat("  Output: TB_NNV_BENCHMARK_OPTIMISTIC.csv\n")
cat("================================================================\n")
