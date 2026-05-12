################################################################################
#
#  TB VACCINATION MODEL — EFFECTIVE COVERAGE OVER TIME (PLAUSIBLE + OPTIMISTIC)
#
################################################################################
#
# PURPOSE:
#   Computes effective vaccine coverage (fraction of Mtb-infected individuals
#   in vaccinated compartments) over time for both scenarios. Three views per
#   scenario:
#     1. Overall:   across all four strata combined
#     2. Stratum:   within each stratum separately
#     3. Targeted:  pooled across targeted strata only (where theta > 0)
#
#   Plausible scenario:  VE = 50%, strategy-specific psi
#   Optimistic scenario: VE = 70%, psi = 50%/yr for all strategies
#
# REQUIRES:
#   - deSolve package
#   - TB_VACCINATION_MODEL_COMPLETE.R in working directory
#
# USAGE:
#   1. Set working directory to folder containing the model script
#   2. Source this file: source("TB_EFFECTIVE_COVERAGE.R")
#   3. Results saved as CSVs in working directory
#
#   To run only one scenario, comment out the unwanted entry in the
#   `scenarios` list below.
#
# OUTPUT FILES:
#   Plausible:
#     - TB_EFFECTIVE_COVERAGE_OVERALL.csv
#     - TB_EFFECTIVE_COVERAGE_STRATUM.csv
#     - TB_EFFECTIVE_COVERAGE_TARGETED.csv
#   Optimistic:
#     - TB_EFFECTIVE_COVERAGE_OVERALL_OPTIMISTIC.csv
#     - TB_EFFECTIVE_COVERAGE_STRATUM_OPTIMISTIC.csv
#     - TB_EFFECTIVE_COVERAGE_TARGETED_OPTIMISTIC.csv
#
# EXPECTED RUNTIME: ~4-6 minutes total
#
# DATE: May 2026
#
################################################################################

suppressPackageStartupMessages({
  library(deSolve)
})

cat("================================================================\n")
cat("  TB VACCINATION MODEL — EFFECTIVE COVERAGE (BOTH SCENARIOS)\n")
cat("================================================================\n\n")

# --- Source the base model (runs calibration automatically) ---
cat("Loading base model and running calibration...\n")
invisible(capture.output(suppressMessages(
  source("TB_VACCINATION_MODEL_COMPLETE.R")
)))
cat("Base model loaded.\n")
cat(sprintf("  Baseline beta = %.4f\n\n", calibration$beta))

baseline_cal <- calibration

# ---- Common simulation parameters ----------------------------------------
sim_years    <- 30
time_points  <- c(1, 2, 3, 5, 10, 15, 20, 25, 30)
stratum_labs <- c("HIV", "Medical", "NUSB", "USB")

# ---- Strategy definitions per scenario -----------------------------------
strategies_plaus <- list(
  all_high_risk = list(name = "All High-Risk (HIV+Med+FB)", short = "All High-Risk",
                       psi = 0.05, theta = c(1, 1, 1, 0)),
  plwh_all_nusb = list(name = "PLWH + All NUSB",            short = "PLWH+NUSB",
                       psi = 0.05, theta = c(1, 0.18, 1, 0)),
  universal     = list(name = "Universal",                  short = "Universal",
                       psi = 0.02, theta = c(1, 1, 1, 1)),
  hiv_medical   = list(name = "HIV + Medical",              short = "HIV+Medical",
                       psi = 0.05, theta = c(1, 1, 0, 0)),
  medical       = list(name = "Medical Only",               short = "Medical",
                       psi = 0.05, theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",                   short = "HIV",
                       psi = 0.10, theta = c(1, 0, 0, 0))
)

OPT_PSI <- 0.50
strategies_opt <- list(
  all_high_risk = list(name = "All High-Risk (HIV+Med+FB)", short = "All High-Risk",
                       psi = OPT_PSI, theta = c(1, 1, 1, 0)),
  plwh_all_nusb = list(name = "PLWH + All NUSB",            short = "PLWH+NUSB",
                       psi = OPT_PSI, theta = c(1, 0.18, 1, 0)),
  universal     = list(name = "Universal",                  short = "Universal",
                       psi = OPT_PSI, theta = c(1, 1, 1, 1)),
  hiv_medical   = list(name = "HIV + Medical",              short = "HIV+Medical",
                       psi = OPT_PSI, theta = c(1, 1, 0, 0)),
  medical       = list(name = "Medical Only",               short = "Medical",
                       psi = OPT_PSI, theta = c(0, 1, 0, 0)),
  hiv           = list(name = "HIV Only",                   short = "HIV",
                       psi = OPT_PSI, theta = c(1, 0, 0, 0))
)

# ---- Define scenarios ----------------------------------------------------
scenarios <- list(
  plausible = list(
    label       = "PLAUSIBLE",
    strategies  = strategies_plaus,
    VE          = VE_d,           # 0.50 from base model
    csv_overall = "TB_EFFECTIVE_COVERAGE_OVERALL.csv",
    csv_stratum = "TB_EFFECTIVE_COVERAGE_STRATUM.csv",
    csv_target  = "TB_EFFECTIVE_COVERAGE_TARGETED.csv"
  ),
  optimistic = list(
    label       = "OPTIMISTIC",
    strategies  = strategies_opt,
    VE          = 0.70,
    csv_overall = "TB_EFFECTIVE_COVERAGE_OVERALL_OPTIMISTIC.csv",
    csv_stratum = "TB_EFFECTIVE_COVERAGE_STRATUM_OPTIMISTIC.csv",
    csv_target  = "TB_EFFECTIVE_COVERAGE_TARGETED_OPTIMISTIC.csv"
  )
)

# ---- Helper functions ----------------------------------------------------
calc_eff_coverage_overall <- function(state) {
  comp     <- extract_compartments(state)
  total_vax  <- sum(comp$L_fv) + sum(comp$L_sv)
  total_ltbi <- sum(comp$L_f) + sum(comp$L_s) + total_vax
  if (total_ltbi > 0) return(total_vax / total_ltbi * 100)
  return(0)
}

calc_eff_coverage_stratum <- function(state, stratum_idx) {
  idx <- (stratum_idx - 1) * 8
  L_f  <- state[idx + 2]
  L_s  <- state[idx + 3]
  L_fv <- state[idx + 4]
  L_sv <- state[idx + 5]
  total <- L_f + L_s + L_fv + L_sv
  vax   <- L_fv + L_sv
  if (total > 0) return(vax / total * 100)
  return(0)
}

calc_eff_coverage_targeted <- function(state, theta_vec) {
  total_ltbi <- 0
  total_vax  <- 0
  for (i in 1:4) {
    if (theta_vec[i] > 0) {
      idx <- (i - 1) * 8
      L_f  <- state[idx + 2]
      L_s  <- state[idx + 3]
      L_fv <- state[idx + 4]
      L_sv <- state[idx + 5]
      total_ltbi <- total_ltbi + L_f + L_s + L_fv + L_sv
      total_vax  <- total_vax + L_fv + L_sv
    }
  }
  if (total_ltbi > 0) return(total_vax / total_ltbi * 100)
  return(0)
}

################################################################################
# RUN EACH SCENARIO
################################################################################

for (scen_name in names(scenarios)) {

  scen <- scenarios[[scen_name]]

  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("#  SCENARIO: %s  (VE = %.0f%%)\n", scen$label, scen$VE * 100))
  cat("################################################################\n\n")

  results_overall  <- data.frame()
  results_stratum  <- data.frame()
  results_targeted <- data.frame()

  for (s in names(scen$strategies)) {
    strat <- scen$strategies[[s]]
    cat(sprintf("  %-15s  (psi=%.0f%%/yr) ... ",
                strat$short, strat$psi * 100))

    params_int <- list(
      beta     = baseline_cal$beta,
      kappa    = c(1, 1, 1, 1),
      psi      = strat$psi,
      theta    = strat$theta,
      mtb_prev = baseline_cal$mtb_prev,
      M        = M_default,
      N        = N_vec,
      VE       = scen$VE
    )

    traj <- run_model(params_int, years = sim_years)

    for (t in time_points) {
      row_idx <- which.min(abs(traj$trajectory[, 1] - t))
      state   <- as.numeric(traj$trajectory[row_idx, -1])

      # Overall
      results_overall <- rbind(results_overall, data.frame(
        year     = t,
        strategy = strat$short,
        eff_coverage_pct = calc_eff_coverage_overall(state)
      ))

      # Stratum-specific
      for (si in 1:4) {
        results_stratum <- rbind(results_stratum, data.frame(
          year     = t,
          strategy = strat$short,
          stratum  = stratum_labs[si],
          eff_coverage_pct = calc_eff_coverage_stratum(state, si),
          targeted = strat$theta[si] > 0
        ))
      }

      # Targeted strata pooled
      results_targeted <- rbind(results_targeted, data.frame(
        year     = t,
        strategy = strat$short,
        eff_coverage_pct = calc_eff_coverage_targeted(state, strat$theta)
      ))
    }

    eq_val <- calc_eff_coverage_overall(traj$final_state)
    cat(sprintf("done (year %d overall: %.1f%%)\n", sim_years, eq_val))
  }

  # Add eff_coverage_plot column for non-targeted masking
  results_stratum$eff_coverage_plot <- ifelse(
    results_stratum$targeted,
    results_stratum$eff_coverage_pct, NA
  )

  # Write CSVs
  write.csv(results_overall,  scen$csv_overall, row.names = FALSE)
  write.csv(results_stratum,  scen$csv_stratum, row.names = FALSE)
  write.csv(results_targeted, scen$csv_target,  row.names = FALSE)

  cat(sprintf("\n  Outputs (%s):\n", scen$label))
  cat(sprintf("    %s\n", scen$csv_overall))
  cat(sprintf("    %s\n", scen$csv_stratum))
  cat(sprintf("    %s\n", scen$csv_target))
}

cat("\n\n================================================================\n")
cat("  ALL SCENARIOS COMPLETE\n")
cat("================================================================\n")
