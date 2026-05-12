################################################################################
#
#  TB VACCINATION MODEL — STRATUM-LEVEL IMPACT (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Extracts per-stratum baseline cases, cases prevented, percent reduction,
#   annual vaccinations, and NNV for both scenarios in a single run. Each
#   scenario writes a separate CSV with identical columns.
#
#   Plausible scenario:   strategy-specific psi, VE = 50% (uses globally
#                         defined `strategies` list from the base model)
#   Optimistic scenario:  psi = 0.50/yr for all strategies, VE = 0.70
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_STRATUM_IMPACT.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   - TB_STRATUM_LEVEL_IMPACT.csv             (plausible)
#   - TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv  (optimistic)
#
# EXPECTED RUNTIME: ~2-4 minutes total
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB VACCINATION MODEL — STRATUM-LEVEL IMPACT (BOTH SCENARIOS)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n", calibration$beta))
cat(sprintf("  Actual 2024 cases = %d\n\n", actual_2024_cases))

stratum_names <- c("HIV", "Medical", "FB", "USB")

# ---- Define optimistic strategies (theta same as plausible, psi/VE override) ---
OPT_VE  <- 0.70
OPT_PSI <- 0.50

strategies_opt <- list(
  "All High-Risk"       = list(theta = c(1, 1, 1, 0),       psi = OPT_PSI, VE = OPT_VE),
  "All NUSB"            = list(theta = c(0.19, 0.18, 1, 0), psi = OPT_PSI, VE = OPT_VE),
  "NUSB Stratum Only"   = list(theta = c(0, 0, 1, 0),       psi = OPT_PSI, VE = OPT_VE),
  "Universal"           = list(theta = c(1, 1, 1, 1),       psi = OPT_PSI, VE = OPT_VE),
  "Medical Only"        = list(theta = c(0, 1, 0, 0),       psi = OPT_PSI, VE = OPT_VE),
  "HIV Only"            = list(theta = c(1, 0, 0, 0),       psi = OPT_PSI, VE = OPT_VE),
  "HIV + Medical"       = list(theta = c(1, 1, 0, 0),       psi = OPT_PSI, VE = OPT_VE),
  "PLWH + All NUSB"     = list(theta = c(1, 0.18, 1, 0),    psi = OPT_PSI, VE = OPT_VE),
  "PLWH + NUSB Stratum" = list(theta = c(1, 0, 1, 0),       psi = OPT_PSI, VE = OPT_VE)
)

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label      = "PLAUSIBLE",
    strategies = strategies,       # from base model script
    csv        = "TB_STRATUM_LEVEL_IMPACT.csv"
  ),
  optimistic = list(
    label      = "OPTIMISTIC",
    strategies = strategies_opt,
    csv        = "TB_STRATUM_LEVEL_IMPACT_OPTIMISTIC.csv"
  )
)

# ---- Baseline (no vaccination) — shared across scenarios -----------------
params_base <- list(
  beta = calibration$beta,
  kappa = c(1, 1, 1, 1),
  psi = 0,
  theta = c(0, 0, 0, 0),
  mtb_prev = calibration$mtb_prev,
  M = M_default,
  N = N_vec,
  VE = 0.50
)

baseline_result <- run_model(params_base, years = 50)
baseline_state  <- baseline_result$final_state
baseline_comp   <- extract_compartments(baseline_state)

# ---- Helper: compute per-stratum row for one strategy --------------------
compute_stratum_rows <- function(strat_name, strat, comp) {
  rows <- data.frame()
  for (i in seq_len(4)) {
    N_i <- comp$S[i] + comp$L_f[i] + comp$L_s[i] + comp$L_fv[i] +
           comp$L_sv[i] + comp$I[i] + comp$I_v[i] + comp$R[i]
    inc_rate      <- ((comp$I[i] + comp$I_v[i]) / N_i) * 100000
    stratum_cases <- inc_rate * N_i / 100000

    N_i_base <- baseline_comp$S[i] + baseline_comp$L_f[i] +
                baseline_comp$L_s[i] + baseline_comp$L_fv[i] +
                baseline_comp$L_sv[i] + baseline_comp$I[i] +
                baseline_comp$I_v[i] + baseline_comp$R[i]
    inc_rate_base          <- ((baseline_comp$I[i] + baseline_comp$I_v[i]) /
                               N_i_base) * 100000
    baseline_stratum_cases <- inc_rate_base * N_i_base / 100000

    stratum_prevented <- baseline_stratum_cases - stratum_cases
    stratum_pct_red   <- ifelse(baseline_stratum_cases > 0,
                                stratum_prevented / baseline_stratum_cases * 100,
                                0)
    stratum_vax <- strat$psi * strat$theta[i] * (comp$L_f[i] + comp$L_s[i])
    stratum_nnv <- ifelse(stratum_prevented > 0,
                          stratum_vax / stratum_prevented, NA)

    rows <- rbind(rows, data.frame(
      Scenario            = strat_name,
      Stratum             = stratum_names[i],
      Baseline_cases      = baseline_stratum_cases,
      With_vaccination    = stratum_cases,
      Cases_prevented     = stratum_prevented,
      Pct_reduction       = stratum_pct_red,
      Annual_vaccinations = stratum_vax,
      NNV                 = stratum_nnv,
      Incidence_per_100k  = inc_rate,
      Baseline_incidence  = inc_rate_base
    ))
  }
  rows
}

build_baseline_rows <- function() {
  rows <- data.frame()
  for (i in seq_len(4)) {
    N_i_base <- baseline_comp$S[i] + baseline_comp$L_f[i] +
                baseline_comp$L_s[i] + baseline_comp$L_fv[i] +
                baseline_comp$L_sv[i] + baseline_comp$I[i] +
                baseline_comp$I_v[i] + baseline_comp$R[i]
    inc_rate_base          <- ((baseline_comp$I[i] + baseline_comp$I_v[i]) /
                               N_i_base) * 100000
    baseline_stratum_cases <- inc_rate_base * N_i_base / 100000
    rows <- rbind(rows, data.frame(
      Scenario            = "Baseline",
      Stratum             = stratum_names[i],
      Baseline_cases      = baseline_stratum_cases,
      With_vaccination    = baseline_stratum_cases,
      Cases_prevented     = 0,
      Pct_reduction       = 0,
      Annual_vaccinations = 0,
      NNV                 = NA,
      Incidence_per_100k  = inc_rate_base,
      Baseline_incidence  = inc_rate_base
    ))
  }
  rows
}

baseline_rows <- build_baseline_rows()

################################################################################
# RUN EACH SCENARIO
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s\n", scen$label))
  cat("################################################################\n\n")

  all_stratum <- data.frame()

  for (strat_name in names(scen$strategies)) {
    strat <- scen$strategies[[strat_name]]
    cat(sprintf("  %-30s  (VE=%.0f%%, psi=%.0f%%/yr) ... ",
                strat_name, strat$VE * 100, strat$psi * 100))

    params <- list(
      beta = calibration$beta,
      kappa = c(1, 1, 1, 1),
      psi = strat$psi,
      theta = strat$theta,
      mtb_prev = calibration$mtb_prev,
      M = M_default,
      N = N_vec,
      VE = strat$VE
    )

    result      <- run_model(params, years = 50)
    final_state <- result$final_state
    comp        <- extract_compartments(final_state)

    all_stratum <- rbind(all_stratum,
                         compute_stratum_rows(strat_name, strat, comp))
    cat("done\n")
  }

  # Append baseline rows
  all_stratum <- rbind(all_stratum, baseline_rows)

  write.csv(all_stratum, scen$csv, row.names = FALSE)
  cat(sprintf("\n  Wrote %d rows to %s\n", nrow(all_stratum), scen$csv))
}

cat("\n\n================================================================\n")
cat("  ALL SCENARIOS COMPLETE\n")
cat("================================================================\n")
cat("\n  Output files:\n")
for (scen_name in names(scenarios)) {
  scen <- scenarios[[scen_name]]
  cat(sprintf("    [%s]  %s\n", scen$label, scen$csv))
}
cat("================================================================\n")
